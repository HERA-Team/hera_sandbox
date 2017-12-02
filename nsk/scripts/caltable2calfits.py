#!/usr/bin/env python2.7
"""
caltable2calfits.py
---------------

convert calibration solutions
from CASA K & G gaincal output
and bandpass output in csv files
into gains in calfits file format

Nicholas Kern
Nov. 2017
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

a = argparse.ArgumentParser(description="Turn CASA calibration solutions in {}.csv files from sky_image.py script into .calfits files")

# Required Parameters
a.add_argument("--fname", type=str, help="output calfits filename.", required=True)
a.add_argument("--uv_file", type=str, help="Path to original miriad uv file of data.", required=True)
# Delay Solution Parameters
a.add_argument("--dly_file", type=str, help="Path to .csv file with antenna delay output from sky_image.py (CASA K cal)")
a.add_argument("--TTonly", default=False, action="store_true", help="only store Tip-Tilt slopes of delay solution.")
a.add_argument("--plot_dlys", default=False, action="store_true", help="plot delay solutions across array.")
# Phase Solution Parameters
a.add_argument("--phs_file", type=str, default=None, help="Path to .csv file with phase output from sky_image.py (CASA G cal)")
a.add_argument("--plot_phs", default=False, action="store_true", help="plot phase solutions across array.")
# Bandpass Solution Parameters
a.add_argument("--bp_file", type=str, help="Path to .csv file with antenna complex bandpass output from sky_image.py (CASA bandpass)")
a.add_argument("--bp_flag_frac", type=float, default=0.7, help="at each freq bin, fraction of antennas flagged needed to broadcast flag to all ants.")
a.add_argument("--bp_pass_flags", default=False, action='store_true', help="propagate all bandpass flags")
a.add_argument("--noBPamp", default=False, action='store_true', help="set BP amplitude solutions to zero.")
a.add_argument("--noBPphase", default=False, action='store_true', help="set BP phase solutions to zero.")
a.add_argument('--bp_medfilt', default=False, action='store_true', help="median filter bandpass solutions")
a.add_argument('--medfilt_kernel', default=7, type=int, help="kernel size (channels) for BP median filter")
a.add_argument('--bp_amp_antavg', default=False, action='store_true', help="average bandpass amplitudes across antennas")
a.add_argument('--bp_TTonly', default=False, action='store_true', help="use only tip-tilt phase mode in bandpass solution")
a.add_argument('--plot_bp', default=False, action='store_true', help="plot final bandpass solutions")
# Misc
a.add_argument("--out_dir", default=None, type=str, help="output directory for calfits file. Default is working directory path")
a.add_argument("--overwrite", default=False, action="store_true", help="overwrite output calfits file if it exists")
a.add_argument("--multiply_gains", default=False, action="store_true", help="change gain_convention from divide to multiply.")
a.add_argument('--silence', default=False, action='store_true', help="silence output to stdout")

def echo(message, mtype=0, verbose=True):
    if verbose:
        if mtype == 0:
            print(message)
        elif mtype == 1:
            print('\n{}\n{}'.format(message, '-'*40))

def caltable2calfits(fname, uv_file, dly_file=None, bp_file=None, out_dir=None, phs_file=None, overwrite=False,
                     TTonly=True, plot_dlys=False, plot_phs=False, gain_convention='multiply', plot_bp=False, noBPphase=False,
                     noBPamp=False, bp_TTonly=False, bp_medfilt=False, medfilt_kernel=5, bp_amp_antavg=False, bp_flag_frac=0.3,
                     bp_pass_flags=False, verbose=True):
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
    if dly_file is not None:
        echo("...processing delays", verbose=verbose)
        # get CASA delays and antennas
        dly_ants, dlys = np.loadtxt(dly_file, delimiter=',', unpack=True)
        dly_ants = dly_ants.astype(np.int).tolist()
        dlys *= 1e-9
        dly_antpos = np.array(map(lambda x: antpos[ants.index(x)], dly_ants))
        # keep only limited information
        if TTonly:
            echo("...keeping only TT component of delays", verbose=verbose)
            A = np.vstack([dly_antpos[:, 0], dly_antpos[:, 1]]).T
            fit = np.linalg.pinv(A.T.dot(A)).dot(A.T).dot(dlys)
            dlys = A.dot(fit)
        # get gains
        dly_gains = np.array([np.exp(2j*np.pi*(freqs-freqs.min())*dlys[dly_ants.index(a)]) \
                                if a in dly_ants else np.ones(Nfreqs, dtype=np.complex) for i, a in enumerate(ants)])
        dly_gains = np.repeat(dly_gains[:, :, np.newaxis], Ntimes, axis=2)[:, :, :, np.newaxis]
        # multiply into gains
        gains *= dly_gains

    # process overall phase if available
    if phs_file is not None:
        echo("...processing gain phase", verbose=verbose)
        # get phase antennas and phases
        phs_ants, phases = np.loadtxt(phs_file, delimiter=',', unpack=True)
        phs_ants = phs_ants.astype(np.int).tolist()
        # construct gains
        phs_gains = np.array([np.exp(1j * phases[phs_ants.index(a)]) if a in phs_ants else np.ones(1, dtype=np.complex) for i, a in enumerate(ants)])
        phs_gains = phs_gains[:, np.newaxis, np.newaxis, np.newaxis]
        # mult into gains
        gains *= phs_gains

    # process bandpass if available
    if bp_file is not None:
        echo("...processing bandpass", verbose=verbose)
        # get bandpass and form complex gains
        bp_data = np.loadtxt(bp_file, delimiter=',', dtype=float)
        bp_freqs = bp_data[:, 0] * 1e6
        bp_Nfreqs = len(bp_freqs)
        bp_data = bp_data[:, 1:]
        bp_header = np.loadtxt(bp_file, delimiter=',', dtype=str, comments='$')[0, 1:]
        bp_header = map(lambda x: x.strip(), bp_header)
        bp_data_dict = dict(zip(bp_header, bp_data.T))
        bp_ants = np.unique(map(lambda x: int(''.join([s for s in x if s.isdigit()])), bp_header))
        bp_gains = np.array(map(lambda a: bp_data_dict[str(a)+'r'] + 1j*bp_data_dict[str(a)+'i'] if a in bp_ants else np.ones(bp_data.shape[0], dtype=np.complex), ants))
        bp_flags = np.array(map(lambda a: bp_data_dict[str(a)+'f'] if a in bp_ants else np.zeros(bp_data.shape[0], dtype=np.bool), ants))
        # apply only flags (for all ants) on data that is flagged more than some fraction
        flag_broadcast = (np.sum(bp_flags, axis=0) / float(bp_flags.shape[0])) > bp_flag_frac
        if bp_pass_flags:
            bp_flags = np.repeat(flag_broadcast.reshape(1, -1).astype(np.bool), Nants, 0)
        else:
            bp_flags = bp_flags + np.repeat(flag_broadcast.reshape(1, -1).astype(np.bool), Nants, 0)
        # configure gains and flags shapes
        bp_gains = bp_gains[:, :, np.newaxis, np.newaxis]
        bp_flags = bp_flags[:, :, np.newaxis, np.newaxis]
        # median filter if desired
        if bp_medfilt:
            echo("...median filtering bandpass", verbose=verbose)
            bp_gains = signal.medfilt(bp_gains.real, kernel_size=(1, medfilt_kernel, 1, 1)) + 1j*signal.medfilt(bp_gains.imag, kernel_size=(1, medfilt_kernel, 1, 1))

        # average amplitude across antenna if desired
        if bp_amp_antavg:
            echo("...averaing bandpass amplitude across antennas", verbose=verbose)
            amp_avg = np.median(np.abs(bp_gains), axis=0).reshape(1, -1, 1, 1)
            bp_gains *= amp_avg / np.abs(bp_gains)

        # take only tip-tilt if desired
        if bp_TTonly:
            echo("...distilling bandpass to overall phase, and tip-tilt phase", verbose=verbose)
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
        bp_freq_select = np.array(map(lambda x: np.argmin(np.abs(freqs-x)), bp_freqs))
        gains[:, bp_freq_select, :, :]  *= bp_gains
        flags[:, bp_freq_select, :, :] += bp_flags

    # check filename
    if fname.split('.')[-1] != "calfits":
        fname += ".calfits"

    # make into calfits
    fname = os.path.join(out_dir, fname)
    echo("...writing to calfits {}".format(fname), verbose=verbose)
    make_calfits(fname, gains, freqs, times, jones, ants, flag_array=flags,
                clobber=overwrite, gain_convention=gain_convention)

    # plot dlys
    if plot_dlys:
        fig, ax = plt.subplots(1, figsize=(8,6))
        ax.grid(True)
        dly_max = np.max(np.abs(dlys*1e9))
        dly_min = -dly_max
        cax = ax.scatter(dly_antpos[:, 0], dly_antpos[:, 1], c=dlys*1e9, s=200, cmap='coolwarm', vmin=dly_min, vmax=dly_max)
        cbar = fig.colorbar(cax)
        cbar.set_label("delay [nanosec]", size=16)
        [ax.text(dly_antpos[i,0]+1, dly_antpos[i,1]+2, str(dly_ants[i]), fontsize=12) for i in range(len(dly_ants))]
        ax.set_xlabel("X [meters]", fontsize=14)
        ax.set_ylabel("Y [meters]", fontsize=14)
        ax.set_title("delay solutions for {}".format(os.path.basename(dly_file)), fontsize=10)
        fig.savefig(dly_file+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()

    # plot phs
    if plot_phs:
        fig, ax = plt.subplots(1, figsize=(8,6))
        ax.grid(True)
        phs_max = np.pi
        phs_min = -np.pi
        cax = ax.scatter(dly_antpos[:, 0], dly_antpos[:, 1], c=phases, s=200, cmap='viridis', vmin=phs_min, vmax=phs_max)
        cbar = fig.colorbar(cax)
        cbar.set_label("phase [radians]", size=16)
        [ax.text(dly_antpos[i,0]+1, dly_antpos[i,1]+2, str(dly_ants[i]), fontsize=12) for i in range(len(dly_ants))]
        ax.set_xlabel("X [meters]", fontsize=14)
        ax.set_ylabel("Y [meters]", fontsize=14)
        ax.set_title("phase solutions for {}".format(os.path.basename(dly_file)), fontsize=10)
        fig.savefig(phs_file+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()

    # plot bandpass
    if plot_bp:
        fig, axes = plt.subplots(2, 1, figsize=(12,6))
        fig.subplots_adjust(hspace=0.3)
        # amplitude
        ax = axes[0]
        ax.grid(True)
        ant_sort = np.argsort(ants)
        ant_sort = ant_sort[map(lambda x: x in bp_ants, np.array(ants)[ant_sort])]
        bp_gains[bp_flags.astype(np.bool)] *= np.nan
        p = ax.plot(np.abs(bp_gains[ant_sort]).squeeze().T, marker='.')
        ax.set_xlabel("channel", fontsize=12)
        ax.set_ylabel("amplitude", fontsize=12)
        ax.set_title("bandpass for {}".format(bp_file), fontsize=10)
        # phase
        ax = axes[1]
        ax.grid(True)
        ax.plot(np.angle(bp_gains[ant_sort]).squeeze().T, marker='.')
        ax.set_xlabel("channel", fontsize=12)
        ax.set_ylabel("phase [radians]", fontsize=12)
        lax = fig.add_axes([0.99, 0.1, 0.05, 0.8])
        lax.axis('off')
        lax.legend(p, sorted(ants), ncol=2)
        fig.savefig(bp_file+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()


if __name__ == "__main__":
    args = a.parse_args()

    gain_convention = 'divide'
    if args.multiply_gains:
        gain_convention = 'multiply'

    caltable2calfits(args.fname, args.uv_file, dly_file=args.dly_file, bp_file=args.bp_file, phs_file=args.phs_file,
                    out_dir=args.out_dir, plot_phs=args.plot_phs, noBPamp=args.noBPamp, noBPphase=args.noBPphase,
                    TTonly=args.TTonly, overwrite=args.overwrite, plot_dlys=args.plot_dlys, gain_convention=gain_convention,
                    bp_medfilt=args.bp_medfilt, medfilt_kernel=args.medfilt_kernel, bp_TTonly=args.bp_TTonly, bp_flag_frac=args.bp_flag_frac,
                    bp_amp_antavg=args.bp_amp_antavg, plot_bp=args.plot_bp, bp_pass_flags=args.bp_pass_flags, verbose=(args.silence is False))

