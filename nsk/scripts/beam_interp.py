#!/usr/bin/env python2.7
"""
beam_interp.py
=============

Interpolate healpix beam
at specified sky coordinates
in topocentric frame
"""
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pyuvdata import UVBeam, UVData
from make_calfits import make_calfits
from scipy import interpolate
import numpy as np
import healpy as hp
import argparse
import os
import sys
import copy

ap = argparse.ArgumentParser(description="interpolate healpix beam in TOPO frame")
ap.add_argument("beamfile", type=str, help="path to a pyuvdata.beamfits file")
ap.add_argument("--uv_file", type=str, default=None, help="path to a miriad file to turn attenuation curves into gains in calfits format")
ap.add_argument("--output_fname", default=None, type=str, help="output .npz filename")
ap.add_argument("--output_dir", default=None, type=str, help="output directory path")
ap.add_argument("--theta", default=0, type=float, help="co-latitude in degrees of desired interpolation point")
ap.add_argument("--phi", default=0, type=float, help="longitude in degrees East of N of desired interpolation point(")
ap.add_argument("--freqs", default=None, type=float, nargs='*', help="frequencies in MHz to interpolate onto")
ap.add_argument("--overwrite", default=False, action='store_true', help="overwrite output file")
ap.add_argument("--plot_atten", default=False, action='store_true', help='plot attenuation factor')


def beam_interp(beamfile, theta, phi, freqs=None, uv_file=None, output_fname=None, output_dir=None,
                overwrite=False, plot_atten=False):
    """
    """

    # overwrite
    if output_dir is None:
        output_dir = os.path.dirname(beamfile)

    if output_fname is None:
        output_fname = beamfile + ".atten.npz"
        output_fname = os.path.join(output_dir, os.path.basename(output_fname))
    else:
        output_fname = os.path.join(output_dir, os.path.basename(output_fname))

    if os.path.exists(output_fname) and overwrite is False:
        print("{} exists, not overwriting".format(output_fname))
        sys.exit(0)

    # load beam model
    uvb = UVBeam()
    uvb.read_beamfits(beamfile)

    # extract beams
    beams = uvb.data_array[0, 0, :, :, :]
    beam_freqs = uvb.freq_array.squeeze()

    # setup beam interpolation
    theta_rad  = theta * np.pi / 180.
    phi_rad = phi * np.pi / 180.

    # iterate over polarization
    beam_interp = []
    for i, p in enumerate(uvb.polarization_array):
        temp_beam = []
        # iterate over frequency
        for j, f in enumerate(beam_freqs):
            temp_beam.append(hp.get_interp_val(beams[i, j, :], theta_rad, phi_rad))

        beam_interp.append(temp_beam)

    beam_interp = np.array(beam_interp)

    # re interpolate onto freqs
    if freqs is not None:
        beam_interp = interpolate.interp1d(beam_freqs, beam_interp, fill_value='extrapolate')(freqs * 1e6)

    else:
        freqs = beam_freqs

    # savefile
    print("saving {}".format(output_fname))
    np.savez(output_fname, freqs=freqs, beam_interp=beam_interp, pols=uvb.polarization_array,
             theta=theta, phi=phi, header="beam_interp.shape = (Npols, Nfreqs)")

    # turn into calfits gains
    if uv_file is not None:
        # load uvfile
        uvd = UVData()
        uvd.read_miriad(uv_file)

        # extract relevant parameters
        uv_freqs = uvd.freq_array.squeeze()
        uv_times = np.unique(uvd.time_array)
        uv_jones = uvb.polarization_array
        uv_antpos, uv_ants = uvd.get_ENU_antpos()
        Nants = len(uv_ants)
        Nfreqs = len(uv_freqs)
        Npols = len(uv_jones)
        Ntimes = len(uv_times)
        calfits_fname = os.path.join(output_dir, '.'.join(os.path.basename(output_fname).split('.')[:-1]) + '.calfits')

        # construct gains
        gains = np.ones((Nants, Nfreqs, Ntimes, Npols), dtype=np.complex)

        # interpolate attenuation curves onto uv_freqs, multiply into gains
        uv_atten = interpolate.interp1d(freqs, beam_interp, fill_value='extrapolate')(uv_freqs)
        gains *= uv_atten.T[np.newaxis, :, np.newaxis, :]

        # make calfits
        print("saving {}".format(calfits_fname))
        make_calfits(calfits_fname, gains, uv_freqs, uv_times, uv_jones, uv_ants, clobber=overwrite, gain_convention='divide')

    # plot attenuation
    if plot_atten:
        fig, ax = plt.subplots(1, 1, figsize=(8,6))

        ax.grid()
        ax.set_xlabel("Frequency [MHz]", fontsize=12)
        ax.set_ylabel("Attenuation Factor", fontsize=12)
        lines = []
        label = []
        for i, p in enumerate(uvb.polarization_array):
            l, = ax.plot(freqs/1e6, beam_interp[i, :])
            lines.append(l)
            label.append("{:d} pol: {:.1f} theta, {:.1f} phi".format(p, theta, phi))
        ax.legend(lines, label, ncol=1)
        figname = ".".join(output_fname.split('.')[:-1]) + ".png"
        print("saving {}".format(figname))
        fig.savefig(figname, dpi=150, bbox_inches='tight', pad=0.05)
        plt.close()

if __name__ == "__main__":
    # parse args
    args = ap.parse_args()

    kwargs = copy.copy(args).__dict__
    kwargs.pop("beamfile")

    print args.freqs
    beam_interp(args.beamfile, **kwargs)
