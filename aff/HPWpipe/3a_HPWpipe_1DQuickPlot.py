
# coding: utf-8

# In[ ]:


# Python Standard Library Packages
import os
import glob
import argparse
import datetime

# Community Developed Packages
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

# HERA Collaboration Packages
import hera_pspec as hp


# In[ ]:


now = datetime.datetime.now()


# In[ ]:


parser = argparse.ArgumentParser()
parser.add_argument(
    '-F',
    '--files',
    help='Designate the hdf5 files to be concatenated in time.',
    nargs='*',
    required=True)
parser.add_argument(
    '-W',
    '--wedge',
    help='Turn wedge folding on',
    action='store_true')
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
#     "-F /lustre/aoc/projects/hera/afortino/H1C_IDR2_1/OCRS/2458098/LSThrs_5.0_to_6.0/*.??.HH.hdf5.*.UVP".split())
# dfiles = sorted(glob.glob(args.files[0]))


# In[ ]:


"""Formatting command line arguments:"""
wedge = args.wedge
if args.savepath is None:
    savepath = os.path.dirname(args.files[0])
else:
    savepath = args.savepath
print 'Saving files to:\n{}'.format(savepath)


# In[ ]:


"""Loading metadata:"""
# This will be deprecated once the UVPspec objects supports adding additional attributes
metadata = np.load(os.path.join(os.path.dirname(dfiles[0]), 'metadata.npz'))

FREQrng = metadata['FREQrng'].tolist()
antpos = metadata['antpos'].tolist()
blp_reds = metadata['blp_reds'].tolist()
blpairs = [(tuple(blpair[0]), tuple(blpair[1])) for blpair in metadata['blpairs'].tolist()]
bls_reds = metadata['bls_reds'].tolist()
baselines = metadata['baselines'].tolist()
norms = metadata['norms'].tolist()

LSTrng = metadata['uvd_extra_keywords'].tolist()['LSTrng'].tolist()
JD = metadata['uvd_extra_keywords'].tolist()['JD']
JDt0 = metadata['uvd_extra_keywords'].tolist()['JDt0']
JDtf = metadata['uvd_extra_keywords'].tolist()['JDtf']
numfiles = metadata['uvd_extra_keywords'].tolist()['numfiles']
ext = metadata['uvd_extra_keywords'].tolist()['ext']
xants = metadata['uvd_extra_keywords'].tolist()['xants'].tolist()


# In[ ]:


"""Creating dictionary for converting between pol integers and pol strings:"""
pol_int_to_str = {1: 'pI', 2: 'pQ', 3: 'pU', 4: 'pV', -5: 'XX', -6: 'YY', -7: 'XY', -8: 'YX'}


# In[ ]:


"""Making plots:"""
# Determine how many rows and columns are needed
if (len(dfiles) <= 4) or (len(dfiles) > 8):
    ncols = len(dfiles)
    nrows = 1
else:
    ncols = 4
    nrows = 2

# Initialize the axes objects
f, axes = plt.subplots(
    ncols=ncols,
    nrows=nrows,
    sharex=True,
    sharey=True,
    figsize=(5 * len(dfiles), 2 * len(dfiles)))
plt.subplots_adjust(wspace=0, hspace=0)

# Plot each file
for dfile, ax in zip(dfiles, axes.flatten()):
    # Load in UVPspec objects
    uvp = hp.UVPSpec()
    uvp.read_hdf5(dfile)
    
    # Average the spectra along redundant baseline groups and time
    uvp.average_spectra(blpair_groups=blp_reds.values(), time_avg=True)
    blp_groups = [bl[0] for bl in blp_reds.values()]
    hp.plot.delay_spectrum(
        uvp,
        blp_groups,
        0,
        uvp.pol_array[0],
        fold=wedge,
        delay=False,
        ax=ax,
        component='abs')
    
    ax.set_title(pol_int_to_str[uvp.pol_array[0]], fontsize=20)
f.suptitle('Frequency Channel Range: {} to {} | Excluded Antennae: {}\n    JD {} from {} to {} | {} Files | LST hour from {} to {}'.format(
    FREQrng[0],
    FREQrng[1],
    str(xants),
    JD,
    JDt0,
    JDtf,
    numfiles,
    LSTrng[0],
    LSTrng[1]))

f.savefig('QuickPlot_at_{mon}:{day}_{hr}:{min}:{sec}.pdf'.format(
    mon=now.month,
    day=now.day,
    hr=now.hour,
    min=now.minute,
    sec=now.second))

