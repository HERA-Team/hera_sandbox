#!/usr/bin/env python2.7
"""
dlys2calfits.py
---------------

convert antenna delays and gains
from CASA K & G gaincal output
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
a.add_argument("--dly_file", type=str, help="Path to .csv file with antenna delay output from sky_image.py (CASA K cal)")
a.add_argument("--uv_file", type=str, help="Path to original miriad uv file of data")
a.add_argument("--phs_file", type=str, default=None, help="Optional path to .csv file with phase output from sky_image.py (CASA G cal)")
a.add_argument("--fname", type=str, help="output calfits filename")
a.add_argument("--TTonly", default=False, action="store_true", help="only store Tip-Tilt slopes of delay solution.")
a.add_argument("--out_dir", default=None, type=str, help="output directory for calfits file. Default is dly_file path")
a.add_argument("--overwrite", default=False, action="store_true", help="overwrite output calfits file if it exists")
a.add_argument("--plot_dlys", default=False, action="store_true", help="plot delay solutions across array.")
a.add_argument("--plot_phs", default=False, action="store_true", help="plot phase solutions across array.")
a.add_argument("--divide_gains", default=False, action="store_true", help="change gain_convention from multiply to divide.")
args = a.parse_args()

gain_convention = 'multiply'
if args.divide_gains:
    gain_convention = 'divide'


def dlys2calfits(dly_file, uv_file, fname, out_dir=None, phs_file=None, overwrite=False,
                TTonly=True, plot_dlys=False, plot_phs=False, gain_convention='multiply'):
    """

    """
    # get out_dir
    if out_dir is None:
        out_dir = os.path.dirname(dly_file)

    # load UVData
    uvd = UVData()
    uvd.read_miriad(uv_file)

    # get ants and antpos
    ants = uvd.antenna_numbers.tolist()
    antpos = uvd.get_ENU_antpos(center=True, pick_data_ants=True)

    # get freqs, times, jones
    freqs = uvd.freq_array.squeeze()
    times = uvd.time_array.reshape(uvd.Ntimes, uvd.Nbls)[:, 0]
    Ntimes = len(times)
    jones = uvd.polarization_array

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
    gains = np.exp(-2j*np.pi*(freqs-freqs.min())*casa_dlys.reshape(-1, 1))
    gains = np.repeat(gains[:, :, np.newaxis], Ntimes, axis=2)[:, :, :, np.newaxis]

    # add phs_file
    if phs_file is not None:
        phs_ants, phases = np.loadtxt(phs_file, delimiter=',', unpack=True)
        phs_ants = phs_ants.astype(np.int)
        phases = phases[map(lambda x: x in casa_ants, phs_ants)]
        phs_gains = np.exp(-1j * phases)
        gains *= phs_gains[:, np.newaxis, np.newaxis, np.newaxis]

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
        cax = ax.scatter(casa_antpos[:, 0], casa_antpos[:, 1], c=casa_dlys*1e9, s=80, cmap='coolwarm', vmin=dly_min, vmax=dly_max)
        fig.colorbar(cax, label="delay [ns]")
        [ax.text(casa_antpos[i,0], casa_antpos[i,1], str(casa_ants[i])) for i in range(len(casa_ants))]
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
        cax = ax.scatter(casa_antpos[:, 0], casa_antpos[:, 1], c=phases, s=80, cmap='viridis', vmin=phs_min, vmax=phs_max)
        fig.colorbar(cax, label="phase [radians]")
        [ax.text(casa_antpos[i,0], casa_antpos[i,1], str(casa_ants[i])) for i in range(len(casa_ants))]
        ax.set_xlabel("X [meters]", fontsize=14)
        ax.set_ylabel("Y [meters]", fontsize=14)
        ax.set_title("phase solutions for {}".format(os.path.basename(dly_file)), fontsize=10)
        fig.savefig(phs_file+'.png', dpi=100, bbox_inches='tight', pad=0.05)
        plt.close()

if __name__ == "__main__":
    dlys2calfits(args.dly_file, args.uv_file, args.fname, phs_file=args.phs_file, out_dir=args.out_dir, plot_phs=args.plot_phs,
                    TTonly=args.TTonly, overwrite=args.overwrite, plot_dlys=args.plot_dlys, gain_convention=gain_convention)

