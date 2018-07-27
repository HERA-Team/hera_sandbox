
# coding: utf-8

# In[ ]:


"""This script will convert N*4 uvf5 files into 4 hdf5 files that each represent N files of one polarization."""
import os
from glob import glob
import argparse
import numpy as np
import astropy.units as u
import astropy.coordinates as aco
from pyuvdata import UVData


# In[ ]:


parser = argparse.ArgumentParser()


# In[ ]:


"""Arguments that determine which files from the directory are grabbed."""
parser.add_argument(
    '-F',
    '--files',
    help='Designate the hdf5 files to be concatenated in time.',
    nargs='*',
    required=True)
parser.add_argument(
    '-P',
    '--pols',
    help='Designate which polarizations to analyze.',
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


# In[ ]:


"""Arguments that determine what data from the files are used."""
parser.add_argument(
    '-L',
    '--LSTrange',
    help='Designate which LST in hours that will be analyzed (e.g.: "6.0 7.0").',
    type=float,
    nargs=2,
    required=True)
parser.add_argument(
    '-R',
    '--FREQrange',
    help='Designate which frequency channels that will be analyzed (e.g.: "580 680").',
    type=int,
    nargs=2,
    required=True)
parser.add_argument(
    '-X',
    '--xants',
    help='Designate which antenna numbers that should be excluded from analysis (e.g.: "0 50 98").',
    type=int,
    nargs='*',
    default=list())


# In[ ]:


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
#     -R 530 730 \
#     -X 0 136 50 2 \
#     -S /lustre/aoc/projects/hera/afortino/H1C_IDR2_1/".split())
# files = sorted(glob(args.files[0]))


# In[ ]:


POLS = sorted(args.pols)
DAY = args.day
EXT = args.ext
LSTrange = args.LSTrange
FREQchans = np.arange(args.FREQrange[0], args.FREQrange[-1] + 1)
XANTS = sorted(args.xants)
SAVEPATH = os.path.join(args.savepath, '{EXT}/{DAY}/FREQchans_{FREQ0}_to_{FREQf}/LSThours_{LST0}_to_{LSTf}'.format(
    DAY=DAY,
    LST0=LSTrange[0],
    LSTf=LSTrange[-1],
    FREQ0=FREQchans[0],
    FREQf=FREQchans[-1],
    EXT=EXT))
os.system('mkdir -p {}'.format(SAVEPATH))
print 'Saving files to:\n{}'.format(SAVEPATH)
print 'Polarizations: {}'.format(POLS)
print 'Extension: {}'.format(EXT)
print 'LST Range: Hour {} to Hour {} '.format(LSTrange[0], LSTrange[-1])
print 'Frequency Channels Range : Channel {} to Channel {}'.format(FREQchans[0], FREQchans[-1])
print 'Excluded Antennae: {}'.format(XANTS)


# In[ ]:


pol_files = {pol: sorted([file for file in files if pol in file]) for pol in POLS}


# In[ ]:


print 'Looping through each polarization...'
for pol, dfiles in pol_files.items():
    print pol
    
    print 'Finding the correct files based on the provided LST range...'
    uvd = UVData()
    files = []
    times = []
    for dfile in dfiles:
        uvd.read_uvh5(dfile, read_data=False)
        LSTrads = np.unique(uvd.lst_array * u.rad)
        LSThours = aco.Angle(LSTrads).hour
        LSTindices = np.where(np.logical_and(LSThours >= LSTrange[0], LSThours <= LSTrange[-1]))[0]

        if LSTindices.size > 0:
            JDtimes = np.take(np.unique(uvd.time_array), LSTindices)
            files.append(dfile)
            times.append(JDtimes.tolist())

    print 'Loading in the correct data based on the provided LST range...'
    uvd = UVData()
    uvd.read_uvh5(
        files[0],
        ant_str='cross',
        freq_chans=FREQchans,
        times=times[0])
    for file, time in zip(files[1:], times[1:]):
        uvdi = UVData()
        uvdi.read_uvh5(
            file, 
            ant_str='cross',
            freq_chans=FREQchans,
            times=time)
        uvd += uvdi

    print 'Making new metadata...'
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

    print 'Saving new metadata...'
    uvd.extra_keywords['JDt0'] = JDt0
    uvd.extra_keywords['JDtf'] = JDtf
    uvd.extra_keywords['JD'] = JD
    uvd.extra_keywords['numfiles'] = numfiles
    uvd.extra_keywords['EXT'] = EXT
    uvd.extra_keywords['FREQchans'] = FREQchans
    uvd.extra_keywords['LSTrange'] = LSTrange
    uvd.extra_keywords['XANTS'] = XANTS
    
    print 'Naming new UVdata object...'
    hdf5 = 'zen.{JD}.{JDt0}_{JDtf}.{pol}.HH.hdf5.{ext}'.format(
        JD=JD,
        JDt0=JDt0,
        JDtf=JDtf,
        pol=pol,
        ext=EXT)
    hdf5 = os.path.join(SAVEPATH, hdf5)
    print 'Writing:'
    print hdf5
    uvd.write_uvh5(hdf5, clobber=True)

