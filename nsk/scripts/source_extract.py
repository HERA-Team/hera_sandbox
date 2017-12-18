#!/usr/bin/env python2.7
"""
source_extract.py
========

get image statistics of FITS file(s)
on a specific source
"""
import astropy.io.fits as fits
from astropy import modeling as mod
import argparse
import os
import sys
import shutil
import glob
import numpy as np
import scipy.stats as stats

a = argparse.ArgumentParser(description="Get FITS image statistics around source at center of image")

a.add_argument("files", type=str, nargs='*', help="filename(s) or glob-parseable string of filename(s)")
a.add_argument("--source", type=str, help="source name, with a <source>.loc file in working directory")
a.add_argument("--radius", type=float, default=2, help="radius in degrees around estimated source position to get source peak")
a.add_argument("--outdir", type=str, default=None, help="output directory")
a.add_argument("--ext", type=str, default=".spectrum.tab", help='extension string for spectrum file')
a.add_argument("--overwrite", default=False, action='store_true', help='overwite output')
a.add_argument("--gaussfit_mult", default=1.0, type=float, help="gaussian fit mask area is gaussfit_mult * synthesized_beam")

def source_extract(imfile, source, radius=1, gaussfit_mult=1.0, **kwargs):

    # open fits file
    hdu = fits.open(imfile)

    # get header and data
    head = hdu[0].header
    data = hdu[0].data.squeeze()

    # calculate beam area in degrees^2
    beam_area = (head["BMAJ"] * head["BMIN"] * np.pi / 4 / np.log(2))

    # calculate pixel area in degrees ^2
    pixel_area = np.abs(head["CDELT1"] * head["CDELT2"])
    Npix_beam = beam_area / pixel_area
    
    # get ra dec coordiantes
    ra_axis = np.linspace(head["CRVAL1"]-head["CDELT1"]*head["NAXIS1"]/2, head["CRVAL1"]+head["CDELT1"]*head["NAXIS1"]/2, head["NAXIS1"])
    dec_axis = np.linspace(head["CRVAL2"]-head["CDELT2"]*head["NAXIS2"]/2, head["CRVAL2"]+head["CDELT2"]*head["NAXIS2"]/2, head["NAXIS2"])
    RA, DEC = np.meshgrid(ra_axis, dec_axis)

    # get source coordinates
    ra, dec = np.loadtxt('{}.loc'.format(source), dtype=str)
    ra, dec = map(float, ra.split(':')), map(float,dec.split(':'))
    ra = (ra[0] + ra[1]/60. + ra[2]/3600.) * 15
    dec = (dec[0] + np.sign(dec[0])*dec[1]/60. + np.sign(dec[0])*dec[2]/3600.)

    # get radius coordinates
    R = np.sqrt((RA - ra)**2 + (DEC - dec)**2)

    # select pixels
    select = R < radius

    # get peak brightness within pixel radius
    peak = np.max(data[select])

    # get rms outside of pixel radius
    rms = np.sqrt(np.mean(data[~select]**2))

    # get peak error
    peak_err = rms / np.sqrt(Npix_beam / 2.0)

    # get frequency of image
    freq = head["CRVAL3"]

    ## fit a 2D gaussian and get integrated and peak flux statistics ##
    # recenter R array by peak flux point and get thata T array
    peak_ind = np.argmax(data[select])
    peak_ra = RA[select][peak_ind]
    peak_dec = DEC[select][peak_ind]
    X = (RA - peak_ra)
    Y = (DEC - peak_dec)
    R = np.sqrt(X**2 + Y**2)
    X[np.where(np.isclose(X, 0.0))] = 1e-5
    T = np.arctan(Y / X)

    # use synthesized beam as data mask
    ecc = head["BMAJ"] / head["BMIN"]
    beam_theta = head["BPA"] * np.pi / 180 + np.pi/2
    EMAJ = R * np.sqrt(np.cos(T+beam_theta)**2 + ecc**2 * np.sin(T+beam_theta)**2)
    fit_mask = EMAJ < (head["BMAJ"] / 2 * gaussfit_mult)
    masked_data = data.copy()
    masked_data[~fit_mask] = 0.0

    # fit 2d gaussian
    gauss_init = mod.functional_models.Gaussian2D(peak, ra, dec, x_stddev=head["BMAJ"]/2, y_stddev=head["BMIN"]/2) 
    fitter = mod.fitting.LevMarLSQFitter()
    gauss_fit = fitter(gauss_init, RA[fit_mask], DEC[fit_mask], data[fit_mask])

    # get gaussian fit properties
    peak_gauss_flux = gauss_fit.amplitude.value
    P = np.array([X, Y]).T
    beam_theta -= np.pi/2
    Prot = P.dot(np.array([[np.cos(beam_theta), -np.sin(beam_theta)], [np.sin(beam_theta), np.cos(beam_theta)]]))
    gauss_cov = np.array([[gauss_fit.x_stddev.value**2, 0], [0, gauss_fit.y_stddev.value**2]])
    model_gauss = stats.multivariate_normal.pdf(Prot, mean=np.array([0,0]), cov=gauss_cov)
    model_gauss *= gauss_fit.amplitude.value / model_gauss.max()
    int_gauss_flux = np.nansum(model_gauss.ravel()) / Npix_beam

    return peak, peak_err, rms, peak_gauss_flux, int_gauss_flux, freq

if __name__ == "__main__":

    # parse args
    args = a.parse_args()

    # sort files
    files = sorted(args.files)

    # get filename
    if args.outdir is None:
        args.outdir = os.path.dirname(os.path.commonprefix(files))
    output_fname = os.path.join(args.outdir, os.path.basename(os.path.splitext(os.path.commonprefix(files))[0] + args.ext))
    if os.path.exists(output_fname) and args.overwrite is False:
        raise IOError("file {} exists, not overwriting".format(output_fname))

    # iterate over files
    peak_flux = []
    peak_flux_err = []
    peak_gauss_flux = []
    int_gauss_flux = []
    freqs = []
    for i, fname in enumerate(files):
        output = source_extract(fname, **vars(args))
        peak_flux.append(output[0])
        peak_flux_err.append(output[2])
        peak_gauss_flux.append(output[3])
        int_gauss_flux.append(output[4])
        freqs.append(output[5])

    freqs = np.array(freqs)
    peak_flux = np.array(peak_flux)
    peak_flux_err = np.array(peak_flux_err)
    peak_gauss_flux = np.array(peak_gauss_flux)
    int_gauss_flux = np.array(int_gauss_flux)

    data = np.vstack([freqs/1e6, peak_flux, peak_flux_err, peak_gauss_flux, int_gauss_flux]).T

    # save spectrum
    print("...saving {}".format(output_fname))
    np.savetxt(output_fname, data, fmt="%8.5f", header="freqs (MHz)\t peak flux (Jy)\t flux err\t peak gaussfit (Jy)\t integ gaussfit (Jy)", delimiter='\t')





