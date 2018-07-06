
# coding: utf-8

# In[ ]:


"""This script will convert N MIRIAD data files into N hdf5 files."""
# Python Standard Library Packages
import os
import glob
import argparse

# Community Developed Packages
import numpy as np

# HERA Collaboration Packages
from pyuvdata import UVData


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument("-d",
    "--day",
    help="Designate the night of IDR2.1 observation to be analyzed.",
    required=True)
parser.add_argument("-e",
    "--ext",
    help="Designate the file extension of the files to be analyzed.",
    required=True)
args = parser.parse_args()
DAY = args.day
EXT = args.ext


# In[ ]:


# DAY = '2458111'
# EXT = 'uvOCRSD'


# In[ ]:


# Path to MIRIAD files
IDR2_1 = '/lustre/aoc/projects/hera/H1C_IDR2/IDR2_1/{}'.format(DAY)

# Retrieves MIRIAD files
FILE_SKELETON = '*{}*{}'.format(DAY, EXT)
FILE_SKELETON = os.path.join(IDR2_1, FILE_SKELETON)
MIRIAD_DFILES = np.array(sorted(glob.glob(FILE_SKELETON)))

# Saves new hdf5 files
SAVE_PATH = '/lustre/aoc/projects/hera/afortino/{}/zen_{}_1time_1pol_HH_{}_hdf5/'.format(DAY, DAY, EXT)
os.system('mkdir -p {}'.format(SAVE_PATH))


# In[ ]:


for dfile in MIRIAD_DFILES:
    print 'Reading: {}'.format(dfile)
    hdf5 = os.path.join(SAVE_PATH, os.path.basename(dfile)) + '.hdf5'
    uvd = UVData()
    uvd.read_miriad(dfile, ant_str='cross')
    uvd.write_uvh5(hdf5, clobber=True)

