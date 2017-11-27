#!/usr/bin/env python2.7
"""
caltable2calfits.py
---------------

convert calibration solutions
from CASA K & G gaincal output
and bandpass output in caltable format
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

a = argparse.ArgumentParser(description="Turn CASA delays.csv files from sky_image.py script into .calfits files")
a.add_argument("--fname", type=str, help="output calfits filename.", required=True)
a.add_argument("--uv_file", type=str, help="Path to original miriad uv file of data.", required=True)
a.add_argument("--bp_file", type=str, help="Path to .csv file with antenna complex bandpass output from sky_image.py (CASA bandpass)")
a.add_argument("--dly_file", type=str, help="Path to .csv file with antenna delay output from sky_image.py (CASA K cal)")
a.add_argument("--uv_file", type=str, help="Path to original miriad uv file of data")
a.add_argument("--phs_file", type=str, default=None, help="Path to .csv file with phase output from sky_image.py (CASA G cal)")
a.add_argument("--out_dir", default=None, type=str, help="output directory for calfits file. Default is dly_file path")
a.add_argument("--TTonly", default=False, action="store_true", help="only store Tip-Tilt slopes of delay solution.")
a.add_argument("--overwrite", default=False, action="store_true", help="overwrite output calfits file if it exists")
a.add_argument("--plot_dlys", default=False, action="store_true", help="plot delay solutions across array.")
a.add_argument("--plot_phs", default=False, action="store_true", help="plot phase solutions across array.")
a.add_argument("--divide_gains", default=False, action="store_true", help="change gain_convention from multiply to divide.")
a.add_argument("--noBPphase", default=False, action='store_true', help="set BP phase solutions to zero.")

args = a.parse_args()

gain_convention = 'multiply'
if args.divide_gains:
    gain_convention = 'divide'

def caltable2calfits(fname, uv_file, dly_file=None, bp_file=None, out_dir=None, phs_file=None, overwrite=False,
                TTonly=True, plot_dlys=False, plot_phs=False, gain_convention='multiply', noBPphase=False):
    """

    """
    # get out_dir
    if out_dir is None:
        out_dir = "./"

    # load UVData
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

    # construct blank gains
    gains = np.ones((Nants, Nfreqs, Ntimes, 1), dtype=np.complex)

    # process delays if available
    if dly_file is not None:
        # get CASA delays and antennas
        casa_ants, casa_dlys = np.loadtxt(dly_file, delimiter=',', unpack=True)
        casa_ants = casa_ants.astype(np.int)
        casa_dlys *= 1e-9
        casa_antpos = np.array(map(lambda x: antpos[ants.index(x)], casa_ants))
        # keep only limited information
        if TTonly:
            A = np.vstack([casa_antpos[:, 0], casa_antpos[:, 1]]).T
            fit = np.linalg.pinv(A.T.dot(A)).dot(A.T).dot(casa_dlys)
            casa_dlys = A.dot(fit)
        # get gains
        dly_gains = np.exp(-2j*np.pi*(freqs-freqs.min())*casa_dlys.reshape(-1, 1))
        dly_gains = np.repeat(gains[:, :, np.newaxis], Ntimes, axis=2)[:, :, :, np.newaxis]
        # multiply into gains
        gains *= dly_gains


    # process overall phase if available
    if phs_file is not None:
        # get phase antennas and phases
        phs_ants, phases = np.loadtxt(phs_file, delimiter=',', unpack=True)
        phs_ants = phs_ants.astype(np.int)
        phases = phases[map(lambda x: x in casa_ants, phs_ants)]
        # construct gains
        phs_gains = np.exp(-1j * phases)[:, np.newaxis, np.newaxis, np.newaxis]
        # multiply into gains
        gains *= phs_gains

    # process bandpass if available
    if bp_file is not None:
        # get bandpass and form complex gains
        bp_data = np.loadtxt(bp_file, delimiter=',', dtype=float)
        bp_freqs = bp_data[:, 0] * 1e6
        bp_data = bp_data[:, 1:]
        bp_header = np.loadtxt(bp_file, delimiter=',', dtype=str, comments='$')[0, 1:]
        bp_header = map(lambda x: x.strip(), bp_header)
        bp_data_dict = dict(zip(bp_header, bp_data.T))
        bp_ants = np.unique(map(lambda x: int(''.join([s for s in x if s.isdigit()])), bp_header))
        bp_gains = np.array(map(lambda a: bp_data_dict[str(a)+'r'] + 1j*bp_data_dict[str(a)+'i'] if a in bp_ants else np.ones(bp_data.shape[0], dtype=np.complex), ants))
        bp_gains = bp_gains[:, :, np.newaxis, np.newaxis]
        # suppress phase if desired
        if noBPphase:
            bp_gains /= np.exp(1j*np.angle(bp_gains))
        # multiply into gains
        bp_freq_select = np.array(map(lambda x: np.argmin(np.abs(freqs-x)), bp_freqs))
        gains[:, bp_freq_select, :, :]  *= bp_gains

    # check filename
    if fname.split('.')[-1] != "calfits":
        fname += ".calfits"

    # make into calfits
    fname = os.path.join(out_dir, fname)
    make_calfits(fname, gains, freqs, times, jones, casa_ants,
                clobber=overwrite, gain_convention=gain_convention)

    # plot dlys
    if plot_dlys:
        fig, ax = plt.subplots(1, figsize=(8,6))
        ax.grid(True)
        dly_max = np.max(np.abs(casa_dlys*1e9))
        dly_min = -dly_max
        cax = ax.scatter(casa_antpos[:, 0], casa_antpos[:, 1], c=casa_dlys*1e9, s=200, cmap='coolwarm', vmin=dly_min, vmax=dly_max)
        cbar = fig.colorbar(cax)
        cbar.set_label("delay [nanosec]", size=16)
        [ax.text(casa_antpos[i,0]+1, casa_antpos[i,1]+2, str(casa_ants[i]), fontsize=12) for i in range(len(casa_ants))]
        ax.set_xlabel("X [meters]", fontsize=14)
        ax.set_ylabel("Y [meters]", fontsize=14)
        ax.set_title("delay solutions for {}".format(os.path.basename(dly_file)), fontsize=10)
        fig.savefig(dly_file+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()

    # plot phs
    if plot_dlys is True and phs_file is not None:
        fig, ax = plt.subplots(1, figsize=(8,6))
        ax.grid(True)
        phs_max = np.pi
        phs_min = -np.pi
        cax = ax.scatter(casa_antpos[:, 0], casa_antpos[:, 1], c=phases, s=200, cmap='viridis', vmin=phs_min, vmax=phs_max)
        cbar = fig.colorbar(cax)
        cbar.set_label("phase [radians]", size=16)
        [ax.text(casa_antpos[i,0]+1, casa_antpos[i,1]+2, str(casa_ants[i]), fontsize=12) for i in range(len(casa_ants))]
        ax.set_xlabel("X [meters]", fontsize=14)
        ax.set_ylabel("Y [meters]", fontsize=14)
        ax.set_title("phase solutions for {}".format(os.path.basename(dly_file)), fontsize=10)
        fig.savefig(phs_file+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()

if __name__ == "__main__":
    caltable2calfits(args.fname, args.uv_file, dly_file=args.dly_file, bp_file=args.bp_file, phs_file=args.phs_file,
                    out_dir=args.out_dir, plot_phs=args.plot_phs, noBPphase=args.noBPphase,
                    TTonly=args.TTonly, overwrite=args.overwrite, plot_dlys=args.plot_dlys, gain_convention=gain_convention)

