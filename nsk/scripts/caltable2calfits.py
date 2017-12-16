#!/usr/bin/env python2.7
"""
caltable2calfits.py
---------------

convert calibration solutions
from CASA K & G gaincal output
and bandpass output in csv files
into gains in calfits file format

Nicholas Kern
Dec. 2017
"""
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pyuvdata import UVCal, UVData
import numpy as np
from make_calfits import make_calfits
import argparse
import os
import scipy.signal as signal
from sklearn import gaussian_process as gp
import copy

a = argparse.ArgumentParser(description="Turn CASA calibration solutions in {}.csv files from sky_image.py script into .calfits files")

# Required Parameters
a.add_argument("--fname", type=str, help="output calfits filename.", required=True)
a.add_argument("--uv_file", type=str, help="Path to original miriad uv file of data.", required=True)
# Delay Solution Parameters
a.add_argument("--dly_files", type=str, default=None, nargs='*', help="Path to .csv file(s) with antenna delay output from sky_image.py (CASA K cal)")
a.add_argument("--TTonly", default=False, action="store_true", help="only store Tip-Tilt slopes of delay solution.")
a.add_argument("--plot_dlys", default=False, action="store_true", help="plot delay solutions across array.")
# Phase Solution Parameters
a.add_argument("--phs_files", type=str, default=None, nargs='*', help="Path to .csv file(s) with phase output from sky_image.py (CASA G cal)")
a.add_argument("--plot_phs", default=False, action="store_true", help="plot phase solutions across array.")
# Amplitude Solution Parameters
a.add_argument("--amp_files", type=str, default=None, nargs='*',  help="Path to .csv file(s) with amplitude output from sky_image.py (CASA G ampcal)")
a.add_argument("--plot_amp", default=False, action='store_true', help='plot amp solution across array.')
a.add_argument("--gain_amp_antavg", default=False, action='store_true', help="average gain amplitudes across antennas")
# Bandpass Solution Parameters
a.add_argument("--bp_files", type=str, default=None, nargs='*', help="Path to .csv file(s) with antenna complex bandpass output from sky_image.py (CASA bandpass)")
a.add_argument("--bp_flag_frac", type=float, default=0.5, help="at each freq bin, fraction of antennas flagged needed to broadcast flag to all ants.")
a.add_argument("--bp_pass_flags", default=False, action='store_true', help="propagate all bandpass flags")
a.add_argument("--noBPamp", default=False, action='store_true', help="set BP amplitude solutions to zero.")
a.add_argument("--noBPphase", default=False, action='store_true', help="set BP phase solutions to zero.")
a.add_argument('--bp_medfilt', default=False, action='store_true', help="median filter bandpass solutions")
a.add_argument('--medfilt_flag', default=False, action='store_true', help='use median filter to flag bad BP solutions in frequency')
a.add_argument('--medfilt_kernel', default=7, type=int, help="kernel size (channels) for BP median filter")
a.add_argument('--bp_amp_antavg', default=False, action='store_true', help="average bandpass amplitudes across antennas")
a.add_argument('--bp_TTonly', default=False, action='store_true', help="use only tip-tilt phase mode in bandpass solution")
a.add_argument('--plot_bp', default=False, action='store_true', help="plot final bandpass solutions")
a.add_argument('--bp_gp_smooth', default=False, action='store_true', help='smooth bandpass w/ gaussian process. Recommended to precede w/ bp_medfilt.')
a.add_argument('--bp_gp_max_dly', default=1000.0, type=float, help="maximum delay in nanosec allowed for gaussian process fit")
a.add_argument('--bp_gp_nrestart', default=1, type=int, help='number of restarts for GP hyperparameter gradient descent.')
a.add_argument('--bp_gp_thin', default=2, type=int, help="thinning factor for freq bins in GP smooth fit")
# Misc
a.add_argument("--out_dir", default=None, type=str, help="output directory for calfits file. Default is working directory path")
a.add_argument("--overwrite", default=False, action="store_true", help="overwrite output calfits file if it exists")
a.add_argument("--multiply_gains", default=False, action="store_true", help="change gain_convention from divide to multiply.")
a.add_argument('--silence', default=False, action='store_true', help="silence output to stdout")

def echo(message, type=0, verbose=True):
    if verbose:
        if type == 0:
            print(message)
        elif type == 1:
            print('\n{}\n{}'.format(message, '-'*40))

def caltable2calfits(fname, uv_file, dly_files=None, amp_files=None, bp_files=None, out_dir=None, phs_files=None, overwrite=False,
                     TTonly=True, plot_dlys=False, plot_phs=False, gain_convention='multiply', plot_bp=False, noBPphase=False,
                     noBPamp=False, bp_TTonly=False, bp_medfilt=False, medfilt_kernel=5, bp_amp_antavg=False, bp_flag_frac=0.3,
                     bp_pass_flags=False, medfilt_flag=False, bp_gp_smooth=False, bp_gp_max_dly=500.0, bp_gp_nrestart=1, bp_gp_thin=2,
                     plot_amp=False, gain_amp_antavg=False, verbose=True):
    """

    """
    # get out_dir
    if out_dir is None:
        out_dir = "./"

    # load UVData
    echo("...loading uv_file", verbose=verbose)
    uvd = UVData()
    uvd.read_miriad(uv_file)

    # get ants and antpos
    antpos, ants = uvd.get_ENU_antpos(center=True, pick_data_ants=True)
    ants = list(ants)

    # get freqs, times, jones
    freqs = uvd.freq_array.squeeze()
    times = uvd.time_array.reshape(uvd.Ntimes, uvd.Nbls)[:, 0]
    Ntimes = len(times)
    Nfreqs = len(freqs)
    Nants = len(ants)
    jones = uvd.polarization_array

    # construct blank gains and flags
    gains = np.ones((Nants, Nfreqs, Ntimes, 1), dtype=np.complex)
    flags = np.zeros((Nants, Nfreqs, Ntimes, 1), dtype=np.float)

    # process delays if available
    if dly_files is not None:
        echo("...processing delays", verbose=verbose, type=1)
        delays = np.zeros_like(ants, dtype=np.float)
        for dly_file in dly_files:
            echo("...processing {}".format(dly_file))
            # get CASA delays and antennas
            dly_ants, dlys = np.loadtxt(dly_file, delimiter=',', unpack=True)
            dly_ants = dly_ants.astype(np.int).tolist()
            dlys *= 1e-9
            # reorder antennas to be consistant w/ ants array
            dlys = np.array(map(lambda a: dlys[dly_ants.index(a)] if a in dly_ants else 0.0, ants))
            # keep only limited information
            if TTonly:
                echo("...keeping only TT component of delays", verbose=verbose)
                A = np.vstack([antpos[:, 0], antpos[:, 1]]).T
                fit = np.linalg.pinv(A.T.dot(A)).dot(A.T).dot(dlys)
                dlys = A.dot(fit)
            # add to delay array
            delays += dlys

        # turn into complex gains
        dly_gains = np.exp(2j*np.pi*(freqs-freqs.min())*delays.reshape(-1, 1))
        dly_gains = np.repeat(dly_gains[:, :, np.newaxis], Ntimes, axis=2)[:, :, :, np.newaxis]

        # multiply into gains
        gains *= dly_gains

    # process overall phase if available
    if phs_files is not None:
        echo("...processing gain phase", verbose=verbose, type=1)
        phases = np.zeros_like(ants, dtype=np.float)
        for phs_file in phs_files:
            echo("...processing {}".format(phs_file))
            # get phase antennas and phases
            phs_ants, phs = np.loadtxt(phs_file, delimiter=',', unpack=True)
            phs_ants = phs_ants.astype(np.int).tolist()
            # reorder to match ants
            phs = np.array([phs[phs_ants.index(a)] if a in phs_ants else 0.0 for a in ants])
            # add to phases
            phases += phs

        # construct gains
        phase_gains = np.exp(1j * phases)
        phase_gains = phase_gains[:, np.newaxis, np.newaxis, np.newaxis]

        # mult into gains
        gains *= phase_gains

    # process overall amplitude if available
    if amp_files is not None:
        echo("...processing gain amp", verbose=verbose, type=1)
        amplitudes = np.ones_like(ants, dtype=np.float)
        for amp_file in amp_files:
            echo("...processing {}".format(amp_file))
            # get amp antenna and amps
            amp_ants, amps = np.loadtxt(amp_file, delimiter=',', unpack=True)
            amp_ants = amp_ants.astype(np.int).tolist()
            # reorder to match ants
            amps = np.array([amps[amp_ants.index(a)] if a in amp_ants else 1.0 for a in ants])
            # add to amplitudes
            amplitudes *= amps

        # average across ants if desired
        if gain_amp_antavg:
            echo("...averaging antenna gain amplitudes", verbose=verbose)
            amplitudes *= np.median(amplitudes) / amplitudes

        # construct gains
        amplitude_gains = amplitudes
        amplitude_gains = amplitude_gains[:, np.newaxis, np.newaxis, np.newaxis]
        # mult into gains
        gains *= amplitude_gains

    # process bandpass if available
    if bp_files is not None:
        echo("...processing bandpass", verbose=verbose, type=1)
        for ii, bp_file in enumerate(bp_files):
            echo("...processing {}".format(bp_file))
            # get bandpass and form complex gains
            bp_data = np.loadtxt(bp_file, delimiter=',', dtype=float)
            bp_freqs = bp_data[:, 0] * 1e6
            bp_Nfreqs = len(bp_freqs)
            bp_data = bp_data[:, 1:]
            bp_header = np.loadtxt(bp_file, delimiter=',', dtype=str, comments='$')[0, 1:]
            bp_header = map(lambda x: x.strip(), bp_header)
            bp_data_dict = dict(zip(bp_header, bp_data.T))
            bp_ants = np.unique(map(lambda x: int(''.join([s for s in x if s.isdigit()])), bp_header))
            bandp_gains = np.array([bp_data_dict[str(a)+'r'] + 1j*bp_data_dict[str(a)+'i'] if a in bp_ants else np.ones(bp_data.shape[0], dtype=np.complex) for a in ants])
            bandp_flags = np.array([bp_data_dict[str(a)+'f'] if a in bp_ants else np.zeros(bp_data.shape[0], dtype=np.bool) for a in ants])
            # apply only flags (for all ants) on data that is flagged more than some fraction
            flag_broadcast = (np.sum(bandp_flags, axis=0) / float(bandp_flags.shape[0])) > bp_flag_frac
            if bp_pass_flags:
                bandp_flags = np.repeat(flag_broadcast.reshape(1, -1).astype(np.bool), Nants, 0)
            else:
                bandp_flags = bandp_flags + np.repeat(flag_broadcast.reshape(1, -1).astype(np.bool), Nants, 0)
            # configure gains and flags shapes
            bandp_gains = bandp_gains[:, :, np.newaxis, np.newaxis]
            bandp_flags = bandp_flags[:, :, np.newaxis, np.newaxis]
            if ii == 0:
                bp_gains = bandp_gains
                bp_flags = bandp_flags
            else:
                bp_gains *= bandp_gains
                bp_flags += bandp_flags

        # median filter if desired
        if bp_medfilt or medfilt_flag:
            echo("...median filtering", verbose=verbose)
            bp_gains_medfilt = signal.medfilt(bp_gains.real, kernel_size=(1, medfilt_kernel, 1, 1)) + 1j*signal.medfilt(bp_gains.imag, kernel_size=(1, medfilt_kernel, 1, 1))

            if medfilt_flag:
                echo("...solving for BP flags w/ medfilt data")
                # get residual and MAD from unfiltered and filtered data
                residual = (np.abs(bp_gains) - np.abs(bp_gains_medfilt)).squeeze()
                resid_std = np.median(np.abs(residual - np.median(residual, axis=1)[:, np.newaxis]), axis=1) * 1.5
                # identify outliers as greater than 10-sigma
                bad = np.array([np.abs(residual[i]) > resid_std[i]*10 for i in range(residual.shape[0])])
                bp_flags += bad[:, :, np.newaxis, np.newaxis]

            if bp_medfilt:
                echo("...bandpass is the median filtered bandpass")
                bp_gains = bp_gains_medfilt

        # average amplitude across antenna if desired
        if bp_amp_antavg:
            echo("...averaging bandpass amplitude across antennas", verbose=verbose)
            amp_avg = np.median(np.abs(bp_gains), axis=0).reshape(1, -1, 1, 1)
            bp_gains *= amp_avg / np.abs(bp_gains)

        # smooth bandpass w/ gaussian process if desired
        if bp_gp_smooth:
            echo("...smoothing with gaussian process", verbose=verbose)
            freq_lambda = 2 * np.pi / (bp_gp_max_dly*1e-9)
            chan_lambda = freq_lambda / np.median(np.diff(freqs)) 
            kernel = 1**2 * gp.kernels.RBF(100.0, (chan_lambda, 2000.0)) + gp.kernels.WhiteKernel(1e-4, (1e-8, 1e0))
            # configure data
            X = np.arange(len(bp_freqs)).reshape(-1, 1)
            bp_gains_real = []
            bp_gains_imag = []
            # iterate over antennas
            for i, a in enumerate(ants):
                if i % 10 == 0: echo("{}/{} ants".format(i, len(ants)), verbose=verbose)
                # skip flagged ants
                if np.isclose(np.nanmean(np.diff(np.abs(bp_gains[i]), axis=0)), 0.0):
                    bp_gains_real.append(np.ones_like(X.squeeze()))
                    bp_gains_imag.append(np.zeros_like(X.squeeze()))
                    continue
                GP = gp.GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=bp_gp_nrestart)
                xdata = X[~bp_flags[i].squeeze()][::bp_gp_thin]
                yreal = bp_gains[i].squeeze()[~bp_flags[i].squeeze()][::bp_gp_thin].real
                yimag = bp_gains[i].squeeze()[~bp_flags[i].squeeze()][::bp_gp_thin].imag
                ydata = np.vstack([yreal, yimag]).T
                GP.fit(xdata, ydata)
                ypred = GP.predict(X)
                bp_gains_real.append(ypred[:, 0])
                bp_gains_imag.append(ypred[:, 1])
            # reconstruct bp gains
            bp_gains = (np.array(bp_gains_real) + 1j*np.array(bp_gains_imag))[:, :, np.newaxis, np.newaxis]

        # take only tip-tilt if desired
        if bp_TTonly:
            echo("...distilling bandpass to only tip-tilt phase", verbose=verbose)
            # get bandpass phase array
            bp_phs = np.angle(bp_gains).reshape(Nants, bp_Nfreqs)
            # form least squares estimate of phase slopes along X and Y
            A = np.vstack([antpos[:, 0], antpos[:, 1]]).T
            fit = np.linalg.pinv(A.T.dot(A)).dot(A.T).dot(bp_phs)
            # median filter
            fit = signal.medfilt(fit, kernel_size=(1, 11))
            # smooth across frequency via FFT
            window = signal.windows.gaussian(fit.shape[1], fit.shape[1]/50.0)
            fit_smooth = signal.fftconvolve(fit, window.reshape(1, -1), mode='same') / np.sum(window)
            # make prediction to get tip-tilt estimates
            bp_TT_phs = A.dot(fit_smooth).reshape(Nants, bp_Nfreqs, 1, 1)
            # update bandpass gains
            bp_gains *= np.exp(1j*bp_TT_phs - 1j*np.angle(bp_gains))

        # suppress amp and/or phase if desired
        if noBPamp:
            echo("...eliminating BP amplitudes", verbose=verbose)
            bp_gains /= np.abs(bp_gains)
        if noBPphase:
            echo("...eliminating BP phases", verbose=verbose)
            bp_gains /= np.exp(1j*np.angle(bp_gains))

        # mult into gains
        bp_freq_select = np.array([np.argmin(np.abs(freqs-x)) for x in bp_freqs])
        gains[:, bp_freq_select, :, :] *= bp_gains
        flags[:, bp_freq_select, :, :] += bp_flags

    # check filename
    if fname.split('.')[-1] != "calfits":
        fname += ".calfits"

    # make into calfits
    fname = os.path.join(out_dir, fname)
    echo("...writing to calfits {}".format(fname), verbose=verbose, type=1)
    make_calfits(fname, gains, freqs, times, jones, ants, flag_array=flags,
                clobber=overwrite, gain_convention=gain_convention)

    # plot dlys
    if plot_dlys:
        fig, ax = plt.subplots(1, figsize=(8,6))
        ax.grid(True)
        dly_max = np.max(np.abs(delays*1e9))
        dly_min = -dly_max
        for i, a in enumerate(ants):
            if a not in dly_ants:
                continue
            cax = ax.scatter(antpos[i, 0], antpos[i, 1], c=delays[i]*1e9, s=200, cmap='coolwarm', vmin=dly_min, vmax=dly_max)
            ax.text(antpos[i,0]+1, antpos[i,1]+2, str(a), fontsize=12)
        cbar = fig.colorbar(cax)
        cbar.set_label("delay [nanosec]", size=16)
        ax.set_xlabel("X [meters]", fontsize=14)
        ax.set_ylabel("Y [meters]", fontsize=14)
        ax.set_title("delay solutions for {}".format(os.path.basename(dly_files[0])), fontsize=10)
        fig.savefig(dly_files[0]+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()

    # plot phs
    if plot_phs:
        fig, ax = plt.subplots(1, figsize=(8,6))
        ax.grid(True)
        phs_max = np.pi
        phs_min = -np.pi
        for i, a in enumerate(ants):
            if a not in phs_ants:
                continue
            cax = ax.scatter(antpos[i, 0], antpos[i, 1], c=phases[i], s=200, cmap='viridis', vmin=phs_min, vmax=phs_max)
            ax.text(antpos[i,0]+1, antpos[i,1]+2, str(a), fontsize=12) 
        cbar = fig.colorbar(cax)
        cbar.set_label("phase [radians]", size=16)
        ax.set_xlabel("X [meters]", fontsize=14)
        ax.set_ylabel("Y [meters]", fontsize=14)
        ax.set_title("phase solutions for {}".format(os.path.basename(phs_files[0])), fontsize=10)
        fig.savefig(phs_files[0]+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()

    # plot amp
    if plot_amp:
        fig, ax = plt.subplots(1, figsize=(8,6))
        ax.grid(True)
        amp_med = np.median(amplitudes[[a in amp_ants for a in ants]])
        amp_std = np.std(amplitudes[[a in amp_ants for a in ants]])
        amp_max = amp_med + amp_std * 2
        amp_min = amp_med - amp_std * 2
        for i, a in enumerate(ants):
            if a not in amp_ants:
                continue
            cax = ax.scatter(antpos[i, 0], antpos[i, 1], c=amplitudes[i], s=200, cmap='rainbow', vmin=amp_min, vmax=amp_max)
            ax.text(antpos[i,0]+1, antpos[i,1]+2, str(a))
        cbar = fig.colorbar(cax)
        cbar.set_label("amplitude", size=16)
        ax.set_xlabel("X [meters]", fontsize=14)
        ax.set_ylabel("Y [meters]", fontsize=14)
        ax.set_title("amplitude solutions for {}".format(os.path.basename(amp_files[0])), fontsize=10)
        fig.savefig(amp_files[0]+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()

    # plot bandpass
    if plot_bp:
        fig, axes = plt.subplots(2, 1, figsize=(12,6))
        fig.subplots_adjust(hspace=0.3)
        # amplitude
        ax = axes[0]
        ax.grid(True)
        bp_gains[bp_flags.astype(np.bool)] *= np.nan
        pls = []
        bp_ant_select = []
        bp_ants = sorted(bp_ants)
        for i, a in enumerate(bp_ants):
            p, = ax.plot(bp_freqs / 1e6, np.abs(bp_gains).squeeze()[ants.index(a)], marker='.')
            pls.append(p)
        ax.set_xlabel("Frequency [MHz]", fontsize=12)
        ax.set_ylabel("Amplitude", fontsize=12)
        ax.set_title("bandpass for {}".format(bp_file), fontsize=10)
        # phase
        ax = axes[1]
        ax.grid(True)
        for i, a in enumerate(bp_ants):
            ax.plot(bp_freqs / 1e6, np.angle(bp_gains).squeeze()[ants.index(a)], marker='.')
        ax.set_xlabel("Frequency [MHz]", fontsize=12)
        ax.set_ylabel("Phase [radians]", fontsize=12)
        lax = fig.add_axes([1.01, 0.1, 0.05, 0.8])
        lax.axis('off')
        ant_sort = np.argsort(ants)
        lax.legend(pls, bp_ants, ncol=2)
        fig.savefig(bp_file+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()


if __name__ == "__main__":
    args = a.parse_args()

    gain_convention = 'divide'
    if args.multiply_gains:
        gain_convention = 'multiply'

    kwargs = copy.copy(args).__dict__
    kwargs.pop('fname')
    kwargs.pop('uv_file')
    kwargs.pop('multiply_gains')
    kwargs['verbose'] = args.silence is False
    kwargs.pop('silence')
    kwargs['gain_convention'] = gain_convention
    caltable2calfits(args.fname, args.uv_file, **kwargs)


