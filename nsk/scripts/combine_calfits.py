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

def combine_calfits(files, fname, outdir=None, overwrite=False, broadcast_flags=True):
    """
    """
    # get io params
    if outdir is None:
        outdir = "./"

    output_fname = os.path.join(outdir, fname)
    if os.path.exists(fname) and overwrite is False:
        raise IOError("{} exists, not overwriting".format(output_fname))

    # iterate over files
    for i, f in enumerate(files):
        if i == 0:
            echo("...loading {}".format(f), verbose=True)
            uvc = UVCal()
            uvc.read_calfits(f)
            f1 = copy.copy(f)

            # set flagged data to unity
            uvc.gain_array[uvc.flag_array] /= uvc.gain_array[uvc.flag_array]

        else:
            uvc2 = UVCal()
            uvc2.read_calfits(f)

            # check matching specs
            if np.isclose(uvc.freq_array, uvc2.freq_array).min() is False:
                print("skipping {} b/c it doesn't match freqs of {}".format(f, f1))
                continue
            elif np.isclose(uvc.jones_array, uvc2.jones_array).min() is False:
                print("skipping {} b/c it doesn't match jones of {}".format(f, f1))
                continue
            elif np.isclose(uvc.time_array, uvc2.time_array).min() is False:
                print("skipping {} b/c it doesn't match times of {}".format(f, f1))
                continue
            elif np.isclose(uvc.spw_array, uvc2.spw_array).min() is False:
                print("skipping {} b/c it doesn't match spw of {}".format(f, f1))
                continue
            elif uvc2.cal_type != uvc.cal_type:
                print("skipping {} b/c its cal_type doesnt match that of {}".format(f, f1))

            # set flagged data to unity
            gain_array = uvc2.gain_array
            gain_array[uvc2.flag_array] /= gain_array[uvc2.flag_array]

            # multiply gain solutions in
            uvc.gain_array *= uvc2.gain_array

            # pass flags
            if broadcast_flags:
                uvc.flag_array += uvc2.flag_array
            else:
                uvc.flag_array = uvc.flag_array * uvc2.flag_array

    # write to file
    echo("...saving {}".format(output_fname))
    uvc.write_calfits(output_fname, clobber=True)


if __name__ == "__main__":
    # parse args
    a = args.parse_args()

    kwargs = dict(vars(a))
    kwargs.pop('files')
    combine_calfits(a.files, **kwargs)

