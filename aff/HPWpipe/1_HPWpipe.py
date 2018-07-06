
# coding: utf-8

# In[ ]:


"""This script will convert N*4 hdf5 files into 4 hdf5 files that each represent N files of one polarization."""
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


# Path to hdf5 files
hdf5_DATA_PATH = '/lustre/aoc/projects/hera/afortino/{}/zen_{}_1time_1pol_HH_{}_hdf5/'.format(DAY, DAY, EXT)

# Retrieves hdf5 files
FILE_SKELETON = '*{}*{}.hdf5'.format(DAY, EXT)
FILE_SKELETON = os.path.join(hdf5_DATA_PATH, FILE_SKELETON)
hdf5_FILES = np.array(sorted(glob.glob(FILE_SKELETON)))

# Saves new hdf5 files
SAVE_PATH = '/lustre/aoc/projects/hera/afortino/{}/zen_{}_times_1pol_HH_{}_hdf5/'.format(DAY, DAY, EXT)
os.system('mkdir -p {}'.format(SAVE_PATH))


# In[ ]:


hdf5_DFILES = {
    'xx': [dfile for dfile in hdf5_FILES if 'xx' in dfile],
    'xy': [dfile for dfile in hdf5_FILES if 'xy' in dfile],
    'yx': [dfile for dfile in hdf5_FILES if 'yx' in dfile],
    'yy': [dfile for dfile in hdf5_FILES if 'yy' in dfile]}


# In[ ]:


SAMPLE_POL = sorted(hdf5_DFILES.keys())[0]

UVD0 = UVData()
UVDf = UVData()
UVD0.read_uvh5(hdf5_DFILES[SAMPLE_POL][0])
UVDf.read_uvh5(hdf5_DFILES[SAMPLE_POL][-1])

JDt0 = np.round(np.unique(UVD0.time_array)[0], 5)
JDtf = np.round(np.unique(UVDf.time_array)[0], 5)
JD = str(int(JDt0))
JDt0 = '{:.5f}'.format(JDt0).split('.')[-1]
JDtf = '{:.5f}'.format(JDtf).split('.')[-1]
numfiles = str(len(hdf5_DFILES[SAMPLE_POL]))

for pol in sorted(hdf5_DFILES):
    print 'Reading: {} files from {}'.format(pol, hdf5_DATA_PATH)
    uvd = UVData()
    uvd.read_uvh5(hdf5_DFILES[pol])
    
    uvd.extra_keywords['JDt0'] = JDt0
    uvd.extra_keywords['JDtf'] = JDtf
    uvd.extra_keywords['JD'] = JD
    uvd.extra_keywords['numfiles'] = numfiles
    uvd.extra_keywords['ext'] = EXT
    
    hdf5 = 'zen.{JD}.{JDt0}_{numfiles}_{JDtf}.{pol}.HH.{ext}.hdf5'.format(
        JD=JD,
        JDt0=JDt0,
        numfiles=numfiles,
        JDtf=JDtf,
        pol=pol,
        ext=EXT)
    hdf5 = os.path.join(SAVE_PATH, hdf5)
    uvd.write_uvh5(hdf5)

