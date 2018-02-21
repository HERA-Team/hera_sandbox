#!/usr/bin/env python2.7
"""
pbcorr.py
=========

Primary Beam Correction
on FITS images, with a 
primary beam healpix model,
or frequency dependent Gaussians.

Nick Kern
December, 2017
"""
import numpy as np
import astropy.io.fits as fits
from astropy import wcs
from pyuvdata import UVBeam
import os
import sys
import glob
import argparse
import shutil
import copy
import healpy
import scipy.stats as stats
import pytz
import datetime
import ephem
from scipy import interpolate


args = argparse.ArgumentParser(description="Primary beam correction on FITS image files, given primary beam model")

args.add_argument("fitsfiles", type=str, nargs='*', help='path of image FITS file(s) to PB correct')

# PB args
args.add_argument("--multiply", default=False, action='store_true', help='multiply data by primary beam, rather than divide')
args.add_argument("--lon", default=21.42830, type=float, help="longitude of observer in degrees east")
args.add_argument("--lat", default=-30.72152, type=float, help="latitude of observer in degrees north")
args.add_argument("--time", type=str, help='time of image observation in UTC "{yr}/{mon}/{day} {hr}:{min}:{sec}" format')

# HEALpix Beam args
args.add_argument("--beamfile", type=str, help="path to primary beam healpix map in pyuvdata.UVBeam format")
args.add_argument("--pol", type=int, default=-5, help="polarization of healpix maps to use for beam models")

# Gaussian Beam args
args.add_argument("--ew_sig", type=float, default=None, nargs='*',
                 help="if no healpix map provided, array of gaussian beam sigmas (per freq) in the east-west direction")
args.add_argument("--ns_sig", type=float, default=None, nargs='*',
                 help="if no healpix map provided, array of gaussian beam sigmas (per freq) in the north-south direction")
args.add_argument("--gauss_freqs", type=float, default=None, nargs='*',
                  help="if no healpix map provided, array of frequencies (Hz) matching length of ew_sig and ns_sig")

# IO args
args.add_argument("--ext", type=str, default="pbcorr", help='extension for output file')
args.add_argument("--outdir", type=str, default=None, help="output directory, default is path to fitsfile")
args.add_argument("--overwrite", default=False, action='store_true', help='overwrite output files')
args.add_argument("--silence", default=False, action='store_true', help='silence output to stdout')
args.add_argument("--spec_cube", default=False, action='store_true', help='assume all fitsfiles are identical except they each have a single but different frequency')

def echo(message, type=0):
    if verbose:
        if type == 0:
            print(message)
        elif type == 1:
            print('\n{}\n{}'.format(message, '-'*40))


if __name__ == "__main__":

    # parse args
    a = args.parse_args()
    verbose = a.silence == False

    # load pb
    if a.beamfile is not None:
        echo("...loading beamfile {}".format(a.beamfile))
        # load beam
        uvb = UVBeam()
        uvb.read_beamfits(a.beamfile)
        # get beam models and beam parameters
        if a.pol not in uvb.polarization_array:
            raise AttributeError("{} not in {} file polarizations {}".format(a.pol, a.beamfile, uvb.polarization_array))
        pol_ind = np.where(uvb.polarization_array == a.pol)[0][0]
        beam_maps = uvb.data_array[0, 0, pol_ind, :, :]
        beam_freqs = uvb.freq_array.squeeze() / 1e6
        Nbeam_freqs = len(beam_freqs)
        beam_nside = healpy.npix2nside(beam_maps.shape[1])

        # construct beam interpolation function
        def beam_interp_func(theta, phi):
            # convert to radians
            theta = copy.copy(theta) * np.pi / 180.0
            phi = copy.copy(phi) * np.pi / 180.0
            shape = theta.shape
            beam_interp = map(lambda m: healpy.get_interp_val(m, theta.ravel(), phi.ravel(), lonlat=False).reshape(shape), beam_maps)
            return np.array(beam_interp)

    # construct pb
    else:
        # construct beam interpolation function
        echo("...constructing beam from Gaussian models")
        if a.ew_sig is None or a.ns_sig is None:
            raise AttributeError("if beamfile is None, then must feed ew_sig and ns_sig")
        def beam_interp_func(theta, phi):
            beam_interp = []
            psi_ew = theta * np.cos(phi)
            psi_ns = theta * np.sin(phi)
            shape = theta.shape
            for i, (ews, nss) in enumerate(zip(a.ew_sig, a.ns_sig)):
                m = stats.multivariate_normal.pdf(np.array([psi_ew.ravel(), psi_ns.ravel()]).T, mean=np.zeros(2),
                                              cov=np.array([[ews**2, 0],[0, nss**2]]))
                beam_interp.append(m.reshape(shape))
            return np.array(beam_interp)
        beam_freqs = np.array(a.gauss_freqs) / 1e6

    # iterate over FITS files
    for i, ffile in enumerate(a.fitsfiles):

        # create output filename
        if a.outdir is None:
            output_dir = os.path.dirname(ffile)
        else:
            output_dir = a.outdir

        output_fname = os.path.basename(ffile)
        output_fname = os.path.splitext(output_fname)
        output_fname = os.path.join(output_dir, output_fname[0] + '.{}'.format(a.ext) + output_fname[1])

        # check for overwrite
        if os.path.exists(output_fname) and a.overwrite is False:
            raise IOError("{} exists, not overwriting".format(output_fname))

        # load hdu
        echo("...loading {}".format(ffile))
        hdu = fits.open(ffile)

        # get header and data
        head = hdu[0].header
        data = hdu[0].data

        # determine if freq precedes stokes in header
        if head['CTYPE3'] == 'FREQ':
            freq_ax = 3
            stok_ax = 4
        else:
            freq_ax = 4
            stok_ax = 3

        # get axes info
        npix1 = head["NAXIS1"]
        npix2 = head["NAXIS2"]
        nstok = head["NAXIS{}".format(stok_ax)]
        nfreq = head["NAXIS{}".format(freq_ax)]

        # get observer coordinates
        observer = ephem.Observer()
        observer.lat = a.lat * np.pi / 180
        observer.lon = a.lon * np.pi / 180
        observer.date = a.time

        # pointing direction
        point_ra, point_dec = np.array(observer.radec_of(0, np.pi/2)) * 180 / np.pi

        # get WCS
        w = wcs.WCS(ffile)

        # get pixel coordinates
        lon_arr, lat_arr = np.meshgrid(np.arange(npix1), np.arange(npix2))
        lon, lat, s, f = w.all_pix2world(lon_arr.ravel(), lat_arr.ravel(), 0, 0, 0)
        lon = lon.reshape(npix2, npix1)
        lat = lat.reshape(npix2, npix1)
        theta = np.sqrt( (lon - point_ra)**2 + (lat - point_dec)**2 )
        phi = np.arctan2((lat-point_dec), (lon-point_ra)) + np.pi

        # get data frequencies
        if freq_ax == 3:
            data_freqs = w.all_pix2world(0, 0, np.arange(nfreq), 0, 0)[2] / 1e6
        else:
            data_freqs = w.all_pix2world(0, 0, 0, np.arange(nfreq), 0)[3] / 1e6
        Ndata_freqs = len(data_freqs)

        if i == 0 or a.spec_cube is False:
            # evaluate primary beam
            echo("...evaluating PB")
            pb = beam_interp_func(theta, phi)

        # interpolate primary beam onto data frequencies
        echo("...interpolating PB")
        pb_shape = (pb.shape[1], pb.shape[2])
        pb_interp = interpolate.interp1d(beam_freqs, pb.reshape(pb.shape[0], -1).T, fill_value='extrapolate')(data_freqs)
        pb_interp = (pb_interp.T).reshape((Ndata_freqs,) + pb_shape)

        # data shape is [naxis4, naxis3, naxis2, naxis1]
        if freq_ax == 3:
            pb_interp = pb_interp[np.newaxis]
        else:
            pb_interp = pb_interp[:, np.newaxis]

        # divide or multiply by primary beam
        if a.multiply is True:
            echo("...multiplying PB into image")
            data_pbcorr = data * pb_interp
        else:
            echo("...dividing PB into image")
            data_pbcorr = data / pb_interp

        echo("...saving {}".format(output_fname))
        fits.writeto(output_fname, data_pbcorr, head, overwrite=True)

        output_pb = '.'.join(output_fname.split('.')[:-2]) + '.pb.' + output_fname.split('.')[-1]
        echo("...saving {}".format(output_pb))
        fits.writeto(output_pb, pb_interp, head, overwrite=True)

