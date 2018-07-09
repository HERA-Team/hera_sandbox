
# coding: utf-8

# In[ ]:


"""This script will convert N*4 hdf5 files into 4 hdf5 files that each represent N files of one polarization."""
import os
from time import time
from glob import glob
import argparse
import numpy as np
from pyuvdata import UVData


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument(
    '-f',
    '--files',
    help='Designate the hdf5 files to be concatenated in time.',
    nargs='*',
    required=True)
parser.add_argument(
    '-d',
    '--day',
    help='Designate which JD the hdf5 files come from.',
    required=True)
parser.add_argument(
    '-e',
    '--ext',
    help='Designate which file extension (i.e., uvOCRS) the designated hdf5 files are.',
    required=True)
parser.add_argument(
    '-s',
    '--savepath',
    help='Designate the path where the new hdf5 files will be saved.',
    required=True)
parser.add_argument(
    '-p',
    '--pols',
    help='Designate which polarizations to analyze.',
    nargs='*',
    required=True)
# args = parser.parse_args("-f /lustre/aoc/projects/hera/afortino/H1C_IDR2_1/2458111/zen_2458111_1time_1pol_HH_uvOCRS_hdf5/*hdf5 -d 2458111 -e uvOCRS -s /lustre/aoc/projects/hera/afortino/H1C_IDR2_1/ -p xx xy yx yy".split())
args = parser.parse_args()


# In[ ]:


files = np.array(sorted(args.files))
day = args.day
ext = args.ext
savepath = os.path.join(args.savepath, '{day}/zen_{day}_times_1pol_HH_{ext}_hdf5'.format(day=day, ext=ext))
pols = np.array(sorted(args.pols))
os.system('mkdir -p {}'.format(savepath))


# In[ ]:


# Uncomment this when running in jupyter notebook
# files = sorted(glob(args.files[0]))


# In[ ]:


pol_files = {pol: sorted([file for file in files if pol in file]) for pol in pols}


# In[ ]:


UVD0 = UVData()
UVDf = UVData()
UVD0.read_uvh5(pol_files[pols[0]][0])
UVDf.read_uvh5(pol_files[pols[0]][-1])

JDt0 = np.round(np.unique(UVD0.time_array)[0], 5)
JDtf = np.round(np.unique(UVDf.time_array)[0], 5)
JD = str(int(JDt0))
JDt0 = '{:.5f}'.format(JDt0).split('.')[-1]
JDtf = '{:.5f}'.format(JDtf).split('.')[-1]
numfiles = len(pol_files[pols[0]])


# In[ ]:


for pol in pols:
    t0 = time()
    print 'Reading: {} files'.format(pol)
    stack = 0
    while len(pol_files[pol]) > 3:
        num_groups = int(np.ceil(len(pol_files[pol]) / 2.))
        for i, file_group in enumerate(np.array_split(pol_files[pol], num_groups)):
            uvd = UVData()
            uvd.read_uvh5(file_group.tolist())

            temp_hdf5 = os.path.join(savepath, '{}.{}.hdf5'.format(stack, str(i).zfill(3)))
            uvd.write_uvh5(temp_hdf5, clobber=True)
        
        if stack != 0:
            for file in new_files:
                os.remove(file)
            new_files = sorted(glob(os.path.join(savepath, '{}.*.hdf5'.format(stack))))
            pol_files[pol] = new_files
        else:
            new_files = sorted(glob(os.path.join(savepath, '{}.*.hdf5'.format(stack))))
            pol_files[pol] = new_files
        
        stack += 1
    
    uvd = UVData()
    uvd.read_uvh5(pol_files[pol])
    
    uvd.extra_keywords['JDt0'] = JDt0
    uvd.extra_keywords['JDtf'] = JDtf
    uvd.extra_keywords['JD'] = JD
    uvd.extra_keywords['numfiles'] = numfiles
    uvd.extra_keywords['ext'] = ext
    
    hdf5 = 'zen.{JD}.{JDt0}_{numfiles}_{JDtf}.{pol}.HH.{ext}.hdf5'.format(
        JD=JD,
        JDt0=JDt0,
        numfiles=numfiles,
        JDtf=JDtf,
        pol=pol,
        ext=ext)
    hdf5 = os.path.join(savepath, hdf5)
    
    print 'Writing: {}'.format(hdf5)
    uvd.write_uvh5(hdf5)

    tf = time()
    print 'That took {} minutes'.format(np.round((tf - t0) / 60., 2))
    
    for file in new_files:
        os.remove(file)
    
    print 

