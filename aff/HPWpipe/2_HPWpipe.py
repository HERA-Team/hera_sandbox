
# coding: utf-8

# In[ ]:


"""This script will convert 4 hdf5 files into N*4 hdf5 files containing 
UVPspec objects where N is the number of redundnant baseline groups."""
# Python Standard Library Packages
import os
# import copy
import glob
import argparse

# Community Developed Packages
import numpy as np

# HERA Collaboration Packages
import hera_pspec as hp
from hera_pspec.data import DATA_PATH
from pyuvdata import UVData
from hera_cal import redcal


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
hdf5_DATA_PATH = '/lustre/aoc/projects/hera/afortino/{}/zen_{}_times_1pol_HH_{}_hdf5/'.format(DAY, DAY, EXT)

# Retrieves hdf5 files
FILE_SKELETON = '*{}*{}.hdf5'.format(DAY, EXT)
FILE_SKELETON = os.path.join(hdf5_DATA_PATH, FILE_SKELETON)
hdf5_FILES = np.array(sorted(glob.glob(FILE_SKELETON)))

# Saves new hdf5 files
SAVE_PATH = '/lustre/aoc/projects/hera/afortino/{}/zen_{}_times_1pol_PS_{}_hdf5/'.format(DAY, DAY, EXT)
os.system('mkdir -p {}'.format(SAVE_PATH))


# In[ ]:


for dfile in hdf5_FILES:
    print 'Reading: {}'.format(dfile)
    
    # Initialize UVData object to read and contain MIRIAD information
    uvd = UVData()
    uvd.read_uvh5(dfile)
    pol = uvd.get_pols()[0].lower()

    # Apply flags
    uvd.data_array *= np.logical_not(uvd.flag_array)

    # Intialize a cosmology and a beam
    cosmo = hp.conversions.Cosmo_Conversions()
    beamfile = os.path.join(DATA_PATH, 'NF_HERA_Beams.beamfits')
    uvb = hp.pspecbeam.PSpecBeamUV(beamfile, cosmo=cosmo)

    # Convert to cosmological units (mK)
    Jy_to_mK = uvb.Jy_to_mK(np.unique(uvd.freq_array), pol="xx")
    uvd.data_array *= Jy_to_mK[None, None, :, None]

    # Shift data and load datasets
    uvd1 = uvd.select(times=np.unique(uvd.time_array)[:-1:2], inplace=False)
    uvd2 = uvd.select(times=np.unique(uvd.time_array)[1::2], inplace=False)
    ds = hp.PSpecData(dsets=[uvd1, uvd2], wgts=[None, None], beam=uvb)

    # Set visibility units
    ds.dsets[0].vis_units = 'mK'
    ds.dsets[1].vis_units = 'mK'
    
    # Phase data (What does this do?)
    ds.rephase_to_dset(0)
    
    """Categorize baselines into physical separation length"""
    # Setup norm binning
    BIN_WIDTH = 0.3
    NORM_BINS = np.arange(0.0, 10000.0, BIN_WIDTH)

    # Retrieve antenna positions in a dictionary
    antpos = {ant: pos for ant, pos in zip(uvd.get_ENU_antpos()[1], uvd.get_ENU_antpos()[0])}

    # Retrieve antenna pairs and bad antennae
    antpairs = uvd.get_antpairs()
    xants1, xants2 = hp.utils.calc_reds(uvd1, uvd2)[3:]
    xants = np.unique(xants1 + xants2)

    # Sort antenna pairs by their physical separation
    reds = {}
    for antpair in antpairs:
        ant0, ant1 = antpair
        if (ant0 in xants) or (ant1 in xants) or (ant0 >= ant1):
            continue
        norm = np.linalg.norm(antpos[ant0] - antpos[ant1])
        norm = np.round(np.digitize(norm, NORM_BINS) * BIN_WIDTH, 1)

        if norm in reds:
            reds[norm].append(antpair)
        else:
            reds[norm] = [antpair]
    norms = sorted(reds.keys())
    
    """Make UVPspec objects for each baseline bin"""
    for norm in norms:
        print 'Making Power Spectrum for baseline group {}m'.format(norm)
        uvp = ds.pspec(
            reds[norm],
            reds[norm],
            (0, 1),
            pols=(pol, pol),
            spw_ranges=SPW,
            taper="blackman-harris",
            verbose=False)
        hdf5 = '{}m.{}.hdf5'.format(norm, pol)
        hdf5 = os.path.join(SAVE_PATH, hdf5)
        uvp.write_hdf5(hdf5, overwrite=True)


# In[ ]:


JDt0 = uvd.extra_keywords['JDt0']
JDtf = uvd.extra_keywords['JDtf']
JD = uvd.extra_keywords['JD']
numfiles = uvd.extra_keywords['numfiles']
EXT = uvd.extra_keywords['ext']

metadata = os.path.join(SAVE_PATH, 'metadata.npz')
np.savez(
    metadata,
    JDt0=JDt0,
    JDtf=JDtf,
    JD=JD,
    numfiles=numfiles,
    EXT=EXT,
    reds=reds,
    norms=norms,
    antpos=antpos,
    xants=xants)

