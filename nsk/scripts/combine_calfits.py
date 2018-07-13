#!/usr/bin/env python2.7
"""
combine_calfits.py
---------------

load multiple calfits gain solutions,
combine them and then write to file.

Nicholas Kern
Dec. 2017
"""
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pyuvdata import UVCal, UVData
import numpy as np
import argparse
import os
import scipy.signal as signal
from sklearn import gaussian_process as gp
import copy
from make_calfits import make_calfits
from hera_cal.utils import combine_calfits

args = argparse.ArgumentParser(description="")

# Required Parameters
args.add_argument("files", type=str, nargs='*', help='path to calfits file(s) to combine (they must be identical in frequency and time structure.')
args.add_argument("--fname", type=str, help="output filename")
# IO Parameters
args.add_argument("--outdir", default=None, type=str, help="output directory")
args.add_argument("--overwrite", default=False, action='store_true', help="overwrite output files")
# Flag Parameters
args.add_argument("--broadcast_flags", default=False, action='store_true', help="broadcast flags of all incoming files to output file")

def echo(message, type=0, verbose=True):
    if verbose:
        if type == 0:
            print(message)
        elif type == 1:
            print('\n{}\n{}'.format(message, '-'*40))

if __name__ == "__main__":
    # parse args
    a = args.parse_args()

    kwargs = dict(vars(a))
    kwargs.pop('files')
    combine_calfits(a.files, **kwargs)

