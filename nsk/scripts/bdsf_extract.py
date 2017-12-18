#!/usr/bin/env python2.7
"""
bdsf_extrac.py
==============

Use bdsf module to fit 2D gaussians and
extract flux properties.

A {source}.loc file contains a single line
following the format

{RA_hour}:{RA_arcmin}:{RA_arcsec}   {Dec_deg}:{Dec_arcmin}:{Dec_arcsec}

where {} are filled in with RA and Dec values
of the source, respectively.

Nick Kern
Nov. 2017
"""
import matplotlib
matplotlib.use('agg')
import bdsf
import os
import sys
import glob
import argparse
import numpy as np
import astropy.io.fits as fits

a = argparse.ArgumentParser(description="Run BSDF source extractor on FITS images.")
a.add_argument("files", type=str, nargs='*', help="space-separated file(s) or glob-parsable filestem of FITS images to run BDSF on.")
a.add_argument("--source", default=None, type=str, help="source name, with {source}.loc file in working directory.")
a.add_argument("--outpath", default=None, type=str, help="directory path of output files, default is path to input file.")
a.add_argument("--extension", default="spectrum.tab", type=str, help="extension of files' common path for final spectrum file of source.")
a.add_argument("--bdsf_extension", default="bdsf.out", type=str, help="extension of filename for bdsf source catalogues")
a.add_argument("--overwrite", default=False, action="store_true", help="overwrite output BDSF catalog file if it exists")
a.add_argument("--trimbox", default=None, type=int, help="number of pixels from image center to search for source")
a.add_argument("--blank_limit", type=float, default=None, help="Jy/beam limit below which pixels are blanked")
a.add_argument("--silence", default=False, action="store_true", help="silence output to stdout")

def echo(message, type=0, verbose=True):
    if verbose:
        if type == 0:
            print(message)
        elif type == 1:
            print('\n{}\n{}'.format(message, '-'*40))

if __name__ == "__main__":
    args = a.parse_args()
    verbose = args.silence is False

    # parse filenames
    filenames = []
    for i, f in enumerate(args.files):
        filenames.extend(glob.glob(f))
    # sort
    filenames = sorted(filenames)
    if len(filenames) == 0:
        raise AttributeError("len(filenames) == 0")

    # run bdsf
    outfiles = []
    for i, fname in enumerate(filenames):
        # get output filename
        if args.outpath is None:
            outpath = os.path.dirname(fname)
        else:
            outpath = args.outpath
        fname_stem = '.'.join(os.path.basename(fname).split('.')[:-1])
        out_fname = os.path.join(outpath, fname_stem + '.' + args.bdsf_extension)

        # check if it exists
        echo("\n::: loading {} :::".format(os.path.basename(fname)))
        if os.path.exists(out_fname) and args.overwrite is False:
            echo("::: {} exists, not overwriting :::".format(out_fname))
            continue

        # configure trimbox
        if args.trimbox is not None:
            with fits.open(fname) as hdu:
                naxis1 = hdu[0].header["NAXIS1"]
                naxis2 = hdu[0].header["NAXIS1"]
                trim_box = (int(naxis1/2)-args.trimbox, int(naxis1/2)+args.trimbox,
                            int(naxis2/2)-args.trimbox, int(naxis2/2)+args.trimbox)
        else:
            trim_box = None
        
        # run bdsf
        try:
            img = bdsf.process_image(fname, quiet=args.silence, trim_box=trim_box, blank_limit=args.blank_limit)
        except RuntimeError:
            continue

        # save catalog to ASCII file
        outfiles.append(out_fname)
        echo("::: saving {} :::".format(out_fname))
        bdsf.output.write_ascii_list(img, filename=out_fname, clobber=args.overwrite)

    # stop if no source provided
    if args.source is None:
        echo("can't proceed with source extraction because args.source is not provided")
        sys.exit(0)

    # get source position in degrees
    ra, dec = np.loadtxt("{}.loc".format(args.source), dtype=str)
    ra, dec = map(float, ra.split(':')), map(float, dec.split(':'))
    ra = (ra[0] + ra[1] / 60. + ra[2] / 3600.) * 180./12.
    dec = (dec[0] + np.sign(dec[0])*dec[1] / 60. + np.sign(dec[0])*dec[2] / 3600.)

    # make filename
    file_stem = os.path.splitext(os.path.basename(os.path.commonprefix(filenames)))[0]
    if args.outpath is None:
        args.outpath = os.path.dirname(os.path.commonprefix(filenames))
    filename = os.path.join(args.outpath, file_stem + '.' + args.extension)
    if os.path.exists(filename) and args.overwrite is False:
        echo("{} exists, not overwriting".format(filename))
        sys.exit(0)

    # setup arrays
    freqs = []
    dist = []
    source_data = []
    # iterate over bdsf catalog files
    for i, fname in enumerate(outfiles):
        if os.path.exists(fname) is False:
            continue
        # open file descriptor
        f = open(fname)
        # read lines
        lines = f.readlines()
        f.close()
        # read values via loadtxt
        data = np.loadtxt(fname, usecols=(4, 6, 8, 9, 10, 11, 16, 18, 20))
        if data.size == 0:
            continue
        elif data.ndim == 1:
            data = data[np.newaxis]
        # append frequency
        freqs.append(float(lines[2].strip().split()[-2]))
        # find the gaussian closest to the source
        ra_dist = data[:, 0] - ra
        de_dist = data[:, 1] - dec
        radius = np.array(map(np.linalg.norm, np.vstack([ra_dist, de_dist]).T))
        s_id = np.argmin(radius)
        # append to arrays
        dist.append(radius[s_id])
        source_data.append(data[s_id])

    freqs = np.array(freqs) / 1e6
    dist = np.array(dist)
    source_data = np.array(source_data)
    output_data = np.concatenate([freqs.reshape(-1, 1), source_data, dist.reshape(-1, 1)], axis=1)

    # write to file
    echo("writing {}".format(filename))

    np.savetxt(filename, output_data, fmt="%10.8f", delimiter='\t',
               header="freq (MHz), RA, Dec, Flux (Jy), Flux Err, Peak Flux (Jy/beam), "\
                      "Peak Flux Err, Maj (deg), Min (deg), PA (Maj deg E of N), source_dist")



