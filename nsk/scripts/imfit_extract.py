"""
imfit_extract.py
===============

Use CASA imfit routine
to fit 2D gaussians to image

"""
import os
import shutil
import glob
import sys
import numpy as np
import argparse

ap = arparse.Argument_Parser(description="Run with casa as: casa -c imfit_extract.py <args>")
ap.add_argument("-c", type=str, help="script name: imfit_extract.py")
ap.add_argument("files", type=str, nargs='*', help="list of files to run imfit over")
ap.add_argument("--region", default=None, type=str, help="path to region file")
ap.add_argument("--filename", default=None, type=str, help="output filename for results")
ap.add_argument("--outdir", default=None, type=str, help="directory path of output files")
ap.add_argument("--excludepix", default="[-1e10, 0]", type=str, help="excludepix paramter of imfit")
ap.add_argument("--overwite", default=False, action='store_true', help="overwrite output files")


if __name__ == "__main__":
    # parse
    a = ap.parse_args()

    # get filename and paths
    files = sorted(files)

    if a.outdir is None:
        outdir = os.path.dirname(os.path.commonprefix(files))

    if a.compfile is None:
        compfile = os.path.basename(os.path.splitext(os.path.commonprefix(files))[0]) + ".comp"

    # iterate over images
    results = []
    for i, im in enumerate(files):
        output = imfit(imagename=im, region=a.region, excludepix=









