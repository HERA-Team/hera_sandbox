
# coding: utf-8

# In[1]:


# Python Standard Library Packages
import os
import copy
import glob

# Community Developed Packages
import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt

# HERA Collaboration Packages
import hera_pspec as hp
from hera_pspec.data import DATA_PATH
from pyuvdata import UVData
from hera_cal import redcal


# In[ ]:


"""Needs to be a step in the pipe to turn into pstokes"""


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
parser.add_argument("-s",
    "--spw",
    help="Designate the spectral window range in channels (i.e., 0_1023).",
    Nargs=*,
    required=True)
args = parser.parse_args()
DAY = args.day
EXT = args.ext
SPW = [(rng.split('_')[0], rng.split('_')[1]) for rng in args.spw]


# In[ ]:


DAY = '2458111'
EXT = 'uvOCRSD'
SPW = [(580, 680)]


# In[ ]:


# Path to hdf5 files
hdf5_DATA_PATH = '/lustre/aoc/projects/hera/afortino/{}/zen_{}_times_1pol_PS_{}_hdf5/'.format(DAY, DAY, EXT)

# Retrieves hdf5 files
FILE_SKELETON = '*{}*{}.hdf5'.format(DAY, EXT)
FILE_SKELETON = os.path.join(hdf5_DATA_PATH, FILE_SKELETON)
hdf5_FILES = np.array(sorted(glob.glob(FILE_SKELETON)))
metadata = os.path.join(hdf5_DATA_PATH, 'metadata.npz')

# Saves new pdf files
SAVE_PATH = '/lustre/aoc/projects/hera/afortino/{}/zen_{}_times_1pol_PS_{}_pdf/'.format(DAY, DAY, EXT)
os.system('mkdir -p {}'.format(SAVE_PATH))


# In[ ]:


metadata = np.load(metadata)
for

