#!/usr/bin/env python2.7
"""
combine_miriad.py
"""
from pyuvdata import UVData
import os
import argparse

args = argparse.ArgumentParser(description="combine multiple miriad files")
args.add_argument("files", type=str, nargs='*', help="space-separated files or glob-parsable search string to concatenate into a single miriad file")
args.add_argument("--fname", type=str, default=None, help="output filename. default is files[0]")
args.add_argument("--outdir", type=str, default=None, help="output directory. default is files[0] directory")
args.add_argument("--overwrite", default=False, action='store_true', help="overwrite even if output file exists")

if __name__ == "__main__":
    a = args.parse_args()

    for i, f in enumerate(a.files):
        uvf = UVData()
        uvf.read_miriad(f)
        if i == 0:
            uvd = uvf
        else:
            uvd += uvf

    if a.fname is None:
        fname = os.path.basename(a.files[0])
        a.overwrite = True
    if a.outdir is None:
        outdir = os.path.dirname(a.files[0])

    out_fname = os.path.join(outdir, fname)
    if os.path.exists(out_fname) and a.overwrite is False:
        print("{} exists, not overwriting".format(out_fname))
    else:
        print("writing {}".format(out_fname))
        uvd.write_miriad(out_fname, clobber=a.overwrite)


