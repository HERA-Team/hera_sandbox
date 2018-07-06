"""This script will convert N MIRIAD data files into N hdf5 files."""
import os
from glob import glob
import argparse
import numpy as np
from pyuvdata import UVData

parser = argparse.ArgumentParser()
parser.add_argument(
    '-f',
    '--files',
    help='Designate the MIRIAD files to be converted, one by one, to hdf5 files.',
    nargs=1,
    required=True)
parser.add_argument(
    '-d',
    '--day',
    help='Designate which JD the MIRIAD files come from.',
    required=True)
parser.add_argument(
    '-e',
    '--ext',
    help='Designate which file extension (i.e., uvOCRS) the designated MIRIAD files are.',
    required=True)
parser.add_argument(
    '-s',
    '--savepath',
    help='Designate the path where the new hdf5 files will be saved.',
    required=True)
args = parser.parse_args()

files = np.array(sorted(glob(args.files)))
day = args.day
ext = args.ext
savepath = os.path.join(args.savepath, '{day}/zen_{day}_1time_1pol_HH_{ext}_hdf5'.format(day=day, ext=ext))
os.system('mkdir -p {}'.format(savepath))

for dfile in files:
    hdf5 = os.path.join(SAVE_PATH, os.path.basename(dfile)) + '.hdf5'
    uvd = UVData()

    print 'Reading: {}'.format(dfile)
    uvd.read_miriad(dfile, ant_str='cross')

    print 'Writing: {}'.format(hdf5)
    uvd.write_uvh5(hdf5, clobber=True)
