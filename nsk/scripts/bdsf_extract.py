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

a = argparse.ArgumentParser(description="Run BSDF source extractor on FITS images.")
a.add_argument("files", type=str, nargs='*', help="space-separated file(s) or glob-parsable filestem of FITS images to run BDSF on.", required=True)
a.add_argument("--source", default=None, type=str, help="source name, with {source}.loc file in working directory.")
a.add_argument("--outpath", default=None, type=str, help="directory path of output files, default is path to input file.")
a.add_argument("--extension", default="spectrum.tab", type=str, help="extension of files' common path for final spectrum file of source.")
a.add_argument("--bdsf_extension", default="bdsf.out", type=str, help="extension of filename for bdsf source catalogues")
a.add_argument("--overwrite", default=False, action="store_true", help="overwrite output BDSF catalog file if it exists")
a.add_argument("--silence", default=False, action="store_true", help="silence output to stdout")

args = a.parse_args()

if __name__ == "__main__":
    # parse filenames
    filenames = []
    for i, f in enumerate(args.files):
        filenames.extend(glob.glob(f))
    # sort
    filenames = sorted(filenames)

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
        outfiles.append(out_fname)

        # check if it exists
        if args.silence is False:
            print("\n::: loading {} :::".format(os.path.basename(fname)))
        if os.path.exists(out_fname) and args.overwrite is False:
            if args.silence is False:
                print("::: {} exists, not overwriting :::".format(out_fname))
            continue

        # run bdsf
        img = bdsf.process_image(fname, quiet=args.silence)

        # save catalog to ASCII file
        if args.silence is False:
            print("::: saving {} :::".format(out_fname))
        bdsf.output.write_ascii_list(img, filename=out_fname, clobber=args.overwrite)

    # stop if no source provided
    if args.source is None:
        print("can't proceed with source extraction because args.source is not provided")
        sys.exit(0)

    # get source position in degrees
    ra, dec = np.loadtxt("{}.loc".format(args.source), dtype=str)
    ra, dec = map(float, ra.split(':')), map(float, dec.split(':'))
    ra = (ra[0] + ra[1] / 60. + ra[2] / 3600.) * 180./12.
    dec = (dec[0] + np.sign(dec[0])*dec[1] / 60. + np.sign(dec[0])*dec[2] / 3600.)

    # make filename
    file_stem = os.path.splitext(os.path.basename(os.path.commonprefix(filenames)))[0]
    filename = os.path.join(args.outpath, file_stem + '.' + args.extension)
    if os.path.exists(filename) and args.overwrite is False:
        print("{} exists, not overwriting".format(filename))
        sys.exit(0)

    # setup arrays
    freqs = []
    dist = []
    source_data = []
    # iterate over bdsf catalog files
    for i, fname in enumerate(out_files):
        if os.path.exists(fname) is False:
            continue
        # open file descriptor
        f = open(fname)
        # read lines
        lines = f.readlines()
        f.close()
        # append frequency
        freqs.append(float(lines[2].strip().split()[-2]))
        # read values via loadtxt
        data = np.loadtxt(fname, usecols=(4, 6, 8, 9, 10, 11, 16, 18, 20))
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
    if args.slience is False:
        print("writing {}".format(filename))

    np.savetxt(filename, output_data, fmt="%10.8f", delimiter='\t',
               header="freq (MHz), RA, Dec, Flux (Jy), Flux Err, Peak Flux (Jy/beam), "\
                      "Peak Flux Err, Maj (deg), Min (deg), PA (Maj deg E of N)")



