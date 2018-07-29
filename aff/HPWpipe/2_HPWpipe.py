
# coding: utf-8

# In[ ]:


"""This script will take one uvf5 file with a UVData object and turn it into multiple
uvf5 files with UVPspec objects in them that each represent the power spectra of a group
of redundant baselines."""
# Python Standard Library Packages
import os
import glob
import argparse

# Community Developed Packages
import numpy as np

# HERA Collaboration Packages
import hera_pspec as hp
from hera_pspec.data import DATA_PATH
from pyuvdata import UVData


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument(
    '-F',
    '--files',
    help='Designate the hdf5 files to be analyzed.',
    nargs='*',
    required=True)
parser.add_argument(
    '-P',
    '--pols',
    help='Designate which polarizations to analyze (e.g.: "pI pQ yx yy xx").',
    nargs='*',
    required=True)
parser.add_argument(
    '-R',
    '--FREQrng',
    help='Designate the frequency range, in channels, that will be analyzed (e.g.: "580 680").',
    type=int,
    nargs=2,
    required=True)
parser.add_argument(
    '-S',
    '--savepath',
    help='Designate the path where the new hdf5 files will be saved. Default is path to data files.')


# In[ ]:


"""Uncomment this code when running as .py:"""
args = parser.parse_args()
dfiles = np.array(sorted(args.files))


# In[ ]:


"""Uncomment this code when running as .ipynb:"""
# args = parser.parse_args(
#     "-F /lustre/aoc/projects/hera/afortino/H1C_IDR2_1/OCRS/{day}/LSThrs_5.0_to_6.0/*OCRS \
#     -P pI pQ pU pV xx xy yx yy \
#     -R 530 730".format(day=day).split())
# dfiles = sorted(glob.glob(args.files[0]))


# In[ ]:


"""Formatting command line arguments:"""
pols = sorted(args.pols)
FREQrng = args.FREQrng
if args.savepath is None:
    savepath = os.path.dirname(args.files[0])
else:
    savepath = args.savepath
savepath = os.path.join(savepath, 'FREQrng_{}_to_{}'.format(FREQrng[0], FREQrng[1]))
os.system('mkdir -p {}'.format(savepath))
print 'Saving files to:\n{}'.format(savepath)
print 'Polarizations: {}'.format(pols)
print 'Frequency Channel Range: {}'.format(FREQrng)


# In[ ]:


"""Defining pol constants:"""
STD_POLS = ['xx', 'xy', 'yx', 'yy']
pS_POLS = ['pI', 'pQ', 'pU', 'pV']


# In[ ]:


"""Loading input files as UVData objects:"""
uvds_std_pols = {std_pol: None for std_pol in STD_POLS}
for dfile in dfiles:
    uvd = UVData()
    uvd.read_uvh5(dfile)
    pol = uvd.get_pols()[0].lower()
    uvds_std_pols[pol] = uvd
    del uvd


# In[ ]:


"""Creating pseudo stokes UVData objects (if requested) and formatting UVdata objects into a dict:"""
uvds = {pol: None for pol in pols}
for pol in pols:

    if pol in STD_POLS:
        uvds[pol] = uvds_std_pols[pol]

    elif pol in pS_POLS:
        if pol == 'pI':
            uvdI = hp.pstokes.construct_pstokes(dset1=uvds_std_pols['xx'], dset2=uvds_std_pols['yy'], pstokes='pI')
            uvds[pol] = uvdI
        if pol == 'pQ':
            uvdQ = hp.pstokes.construct_pstokes(dset1=uvds_std_pols['xx'], dset2=uvds_std_pols['yy'], pstokes='pQ')
            uvds[pol] = uvdQ
        if pol == 'pU':
            uvdU = hp.pstokes.construct_pstokes(dset1=uvds_std_pols['xy'], dset2=uvds_std_pols['yx'], pstokes='pU')
            uvds[pol] = uvdU
        if pol == 'pV':
            uvdV = hp.pstokes.construct_pstokes(dset1=uvds_std_pols['xy'], dset2=uvds_std_pols['yx'], pstokes='pV')
            uvds[pol] = uvdV


# In[ ]:


"""Making UVPspec objects:"""
for pol, uvd in uvds.items():
    # Apply flags
    uvd.data_array *= np.logical_not(uvd.flag_array)

    # Intialize a cosmology and a beam
    if pol in STD_POLS:
        beamfile = os.path.join(DATA_PATH, 'HERA_NF_dipole_power.beamfits')
    elif pol in pS_POLS:
        beamfile = os.path.join(DATA_PATH, 'HERA_NF_pstokes_power.beamfits')
    cosmo = hp.conversions.Cosmo_Conversions()
    uvb = hp.pspecbeam.PSpecBeamUV(beamfile, cosmo=cosmo)

    # Convert to cosmological units (mK)
    if ('C' in uvd.extra_keywords['ext']) or ('K' in uvd.extra_keywords['ext']):
        Jy_to_mK = uvb.Jy_to_mK(np.unique(uvd.freq_array), pol=pol)
        uvd.data_array *= Jy_to_mK[None, None, :, None]

    # Shift data and load datasets
    uvd1 = uvd.select(times=np.unique(uvd.time_array)[:-1:2], inplace=False)
    uvd2 = uvd.select(times=np.unique(uvd.time_array)[1::2], inplace=False)
    ds = hp.PSpecData(dsets=[uvd1, uvd2], wgts=[None, None], beam=uvb)

    # Set visibility units
    if ('C' in uvd.extra_keywords['ext']) or ('K' in uvd.extra_keywords['ext']):
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
    xants = uvd.extra_keywords['xants']

    # Sort antenna pairs by their physical separation
    blpairs, blp_reds = [], {}
    baselines, bls_reds = [], {}
    norms = []
    for antpair in antpairs:
        ant0, ant1 = antpair
        if (ant0 in xants) or (ant1 in xants) or (ant0 >= ant1):
            continue
        baselines.append(antpair)
        blpair = (antpair, antpair)
        blpairs.append(blpair)
        norm = np.linalg.norm(antpos[ant0] - antpos[ant1])
        norm = np.round(np.digitize(norm, NORM_BINS) * BIN_WIDTH, 1)
        norms.append(norm)

        if norm in bls_reds:
            bls_reds[norm].append(antpair)
            blp_reds[norm].append(blpair)
        else:
            bls_reds[norm] = [antpair]
            blp_reds[norm] = [blpair]
    norms = sorted(np.unique(norms))

    """Make UVPspec object"""
    uvp = ds.pspec(
        baselines,
        baselines,
        (0, 1),
        pols=[(pol, pol)],
        spw_ranges=[(FREQrng[0], FREQrng[-1])],
        taper="blackman-harris",
        verbose=False)

    """Name and save UVPspec object"""
    hdf5 = 'zen.{JD}.{JDt0}_{JDtf}.{pol}.HH.hdf5.{ext}.UVP'.format(
        JD=uvd.extra_keywords['JD'],
        JDt0=uvd.extra_keywords['JDt0'],
        JDtf=uvd.extra_keywords['JDtf'],
        pol=pol,
        ext=uvd.extra_keywords['ext'])
    hdf5 = os.path.join(savepath, hdf5)
    print 'Writing:'
    print hdf5
    uvp.write_hdf5(hdf5, overwrite=True)


# In[ ]:


"""Saving additional metadata:"""
metadata = os.path.join(savepath, 'metadata.npz')
np.savez(
    metadata,
    uvd_extra_keywords=uvd.extra_keywords,
    bls_reds=bls_reds,
    blp_reds=blp_reds,
    baselines=baselines,
    blpairs=blpairs,
    norms=norms,
    antpos=antpos,
    FREQrng=FREQrng)

