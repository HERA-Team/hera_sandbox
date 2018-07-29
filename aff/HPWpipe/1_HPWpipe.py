
# coding: utf-8

# In[ ]:


"""This script will concatenate uvf5 files in time based on certain parameters."""
# Python Standard Library Packages
import os
import glob
import argparse

# Community Developed Packages
import numpy as np
import astropy.units as u
import astropy.coordinates as aco

# HERA Collaboration Packages
from pyuvdata import UVData


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument(
    '-F',
    '--files',
    help='Designate the hdf5 files to be concatenated in time.',
    nargs='*',
    required=True)
parser.add_argument(
    '-P',
    '--pols',
    help='Designate which polarizations to analyze (e.g.: "xx xy yx yy").',
    nargs='*',
    required=True)
parser.add_argument(
    '-D',
    '--day',
    help='Designate which JD the hdf5 files come from.',
    required=True)
parser.add_argument(
    '-E',
    '--ext',
    help='Designate which file extension (i.e., uvOCRS) the designated hdf5 files are.',
    required=True)
parser.add_argument(
    '-L',
    '--LSTrng',
    help='Designate which LST in hours that will be analyzed (e.g.: "6.0 7.0").',
    type=float,
    nargs=2,
    required=True)
parser.add_argument(
    '-X',
    '--xants',
    help='Designate which antenna numbers that should be excluded from analysis (e.g.: "0 50 98").',
    type=int,
    nargs='*',
    default=list())
parser.add_argument(
    '-S',
    '--savepath',
    help='Designate the path where the new hdf5 files will be saved.',
    required=True)


# In[ ]:


"""Uncomment this code when running as .py:"""
args = parser.parse_args()
files = np.array(sorted(args.files))


# In[ ]:


"""Uncomment this code when running as .ipynb:"""
# args = parser.parse_args(
#     "-F /lustre/aoc/projects/hera/H1C_IDR2/IDR2_1/2458098/zen.2458098.?????.??.HH.uvh5.OCRS \
#     -P xx xy yx yy\
#     -D 2458098 \
#     -E OCRS \
#     -L 5.0 6.0 \
#     -X 0 136 50 2 \
#     -S /lustre/aoc/projects/hera/afortino/H1C_IDR2_1/".split())
# files = sorted(glob.glob(args.files[0]))


# In[ ]:


"""Formatting command line arguments:"""
pols = sorted(args.pols)
day = args.day
ext = args.ext
LSTrng = args.LSTrng
xants = sorted(args.xants)
savepath = os.path.join(args.savepath, '{ext}/{day}/LSThrs_{LST0}_to_{LSTf}'.format(
    ext=ext,
    day=day,
    LST0=LSTrng[0],
    LSTf=LSTrng[-1]))
os.system('mkdir -p {}'.format(savepath))
print 'Saving files to:\n{}'.format(savepath)
print 'Polarizations: {}'.format(pols)
print 'Extension: {}'.format(ext)
print 'LST Range: Hour {} to Hour {} '.format(LSTrng[0], LSTrng[-1])
print 'Excluded Antennae: {}'.format(xants)


# In[ ]:


"""Looping through each polarization, finding correct files based on given LST range:"""
pol_files = {pol: sorted([file for file in files if pol in file]) for pol in pols}
for pol, dfiles in pol_files.items():
    print pol
    
    # Finding the correct files based on the provided LST range
    uvd = UVData()
    files = []
    times = []
    for dfile in dfiles:
        uvd.read_uvh5(dfile, read_data=False)
        LSTrads = np.unique(uvd.lst_array * u.rad)
        LSThrs = aco.Angle(LSTrads).hour
        LSTindices = np.where(np.logical_and(LSThrs >= LSTrng[0], LSThrs <= LSTrng[-1]))[0]

        if LSTindices.size > 0:
            JDtimes = np.take(np.unique(uvd.time_array), LSTindices)
            files.append(dfile)
            times.append(JDtimes.tolist())

    # Loading in the correct data based on the provided LST range
    uvd = UVData()
    uvd.read_uvh5(
        files[0],
        ant_str='cross',
        times=times[0])
    for file, time in zip(files[1:], times[1:]):
        uvdi = UVData()
        uvdi.read_uvh5(
            file, 
            ant_str='cross',
            times=time)
        uvd += uvdi

    # Making new metadata
    UVD0 = UVData()
    UVDf = UVData()
    UVD0.read_uvh5(files[0], read_data=False)
    UVDf.read_uvh5(files[-1], read_data=False)

    JDt0 = np.round(np.unique(UVD0.time_array)[0], 5)
    JDtf = np.round(np.unique(UVDf.time_array)[0], 5)
    JD = str(int(JDt0))
    JDt0 = '{:.5f}'.format(JDt0).split('.')[-1]
    JDtf = '{:.5f}'.format(JDtf).split('.')[-1]
    numfiles = len(files)

    # Saving new metadata
    uvd.extra_keywords['JDt0'] = JDt0
    uvd.extra_keywords['JDtf'] = JDtf
    uvd.extra_keywords['JD'] = JD
    uvd.extra_keywords['numfiles'] = numfiles
    uvd.extra_keywords['ext'] = ext
    uvd.extra_keywords['LSTrng'] = LSTrng
    uvd.extra_keywords['xants'] = xants
    
    # Naming new UVdata object
    hdf5 = 'zen.{JD}.{JDt0}_{JDtf}.{pol}.HH.hdf5.{ext}'.format(
        JD=JD,
        JDt0=JDt0,
        JDtf=JDtf,
        pol=pol,
        ext=ext)
    hdf5 = os.path.join(savepath, hdf5)
    print 'Writing:'
    print hdf5
    uvd.write_uvh5(hdf5, clobber=True)

