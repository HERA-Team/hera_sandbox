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


args = argparse.ArgumentParser(description="")
args.add_argument("fitsfiles", type=str, nargs='*', help='path to image FITS file(s) to PB correct')
# Beam HEALpix args
args.add_argument("--beamfile", type=str, help="path to primary beam healpix map in pyuvdata.UVBeam format")
args.add_argument("--pol", type=int, default=-5, help="polarization of healpix maps to use for beam models")
# Gaussian Beam args
args.add_argument("--ew_sig", type=float, default=None, nargs='*',
                 help="if no healpix map provided, array of gaussian beam sigma (per freq) in the east-west direction")
args.add_argument("--ns_sig", type=float, default=None, nargs='*',
                 help="if no healpix map provided, array of gaussian beam sigma (per freq) in the north-south direction")
# IO args
args.add_argument("--extension", type=str, default="pbcorr", help='extension for output file')
args.add_argument("--outdir", type=str, default=None, help="output directory, default is path to fitsfile")
args.add_argument("--overwrite", default=False, action='store_true', help='overwrite output files')
args.add_argument("--silence", default=False, action='store_true', help='silence output to stdout')
# Misc args
args.add_argument("--lon", default=21.42830, type=float, help="longitude of observer in degrees east")
args.add_argument("--lat", default=-30.72152, type=float, help="latitude of observer in degrees north")
args.add_argument("--time", type=str, help="time of image observation in UTC {yr}/{mon}/{day} {hr}:{min}:{sec} format")

def echo(message, type=0):
    if verbose:
        if type == 0:
            print(message)
        elif type == 1:
            print('\n{}\n{}'.format(message, '-'*40))


if __name__ == "__main__":

    # parse args
    a = args.parse_args()
    verbose = a.silence is False

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
            theta *= np.pi / 180.0
            phi *= np.pi / 180.0
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

    # iterate over FITS files
    for i, ffile in enumerate(a.fitsfiles):
        echo("...loading {}".format(ffile))

        # load hdu
        hdu = fits.open(ffile)

        # get header and data
        head = hdu[0].header
        data = hdu[0].data
        npix1 = head["NAXIS1"]
        npix2 = head["NAXIS2"]
        nfreq = head["NAXIS3"]
        nstok = head["NAXIS4"]
        obsra = head["OBSRA"]
        obsde = head["OBSDEC"]

        # get observer coordinates
        observer = ephem.Observer()
        observer.lat = a.lat * np.pi / 180
        observer.lon = a.lon * np.pi / 180
        observer.date = a.time
        point_ra, point_dec = np.array(observer.radec_of(0, np.pi/2)) * 180 / np.pi

        # get WCS
        w = wcs.WCS(ffile)

        # get pixel coordinates
        lon, lat, imfreqs, stokes = w.all_pix2world(np.arange(npix1), np.arange(npix2), np.arange(nfreq), np.arange(nstok), 0)








        echo("...saving {}".format(output_fname))



