
# coding: utf-8

# In[2]:


import os
import glob
import argparse
import copy

import numpy as np
import astropy.units as u
import astropy.constants as c
import astropy.coordinates as aco
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import hera_pspec as hp
from hera_pspec.data import DATA_PATH
from pyuvdata import UVData


# In[3]:


parser = argparse.ArgumentParser()

parser.add_argument(
    '-F',
    '--files',
    help='Designate the hdf5 files to be concatenated in time.',
    nargs='*',
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
    '-R',
    '--FREQrng',
    help='Designate the frequency range, in channels, that will be analyzed (e.g.: "580 680").',
    type=int,
    nargs=2,
    required=True)
parser.add_argument(
    '-S',
    '--savepath',
    help='Designate the path where the new hdf5 files will be saved.',
    default='./')

"""Uncomment this code when running as .py"""
# args = parser.parse_args()
# dfiles = np.array(sorted(args.files))

"""Uncomment this code when running as .ipynb"""
args = parser.parse_args(
    "-F /lustre/aoc/projects/hera/H1C_IDR2/IDR2_1/2458098/zen.2458098.?????.xx.HH.uvh5.OCRS \
    -L 5.0 6.0 \
    -X 0 136 50 2 98 137 11 \
    -R 530 730".split())
dfiles = sorted(glob.glob(args.files[0]))

LSTrng = args.LSTrng
xants = sorted(args.xants)
FREQrng = args.FREQrng
savepath = args.savepath
os.system('mkdir -p {}'.format(savepath))
print 'Saving files to:\n{}'.format(savepath)
print 'LST Range: Hour {} to Hour {} '.format(LSTrng[0], LSTrng[-1])
print 'Excluded Antennae: {}'.format(xants)
print 'Frequency Channel Range: {}'.format(FREQrng)


# In[4]:


get_ipython().run_cell_magic(u'time', u'', u'# Finding the correct files based on the provided LST range\nuvd = UVData()\nfiles = []\ntimes = []\nfor dfile in dfiles:\n    uvd.read_uvh5(dfile, read_data=False)\n    LSTrads = np.unique(uvd.lst_array * u.rad)\n    LSThrs = aco.Angle(LSTrads).hour\n    LSTindices = np.where(np.logical_and(LSThrs >= LSTrng[0], LSThrs <= LSTrng[-1]))[0]\n\n    if LSTindices.size > 0:\n        JDtimes = np.take(np.unique(uvd.time_array), LSTindices)\n        files.append(dfile)\n        times.append(JDtimes.tolist())')


# In[5]:


get_ipython().run_cell_magic(u'time', u'', u"# Loading in the correct data based on the provided LST range\nuvd = UVData()\nuvd.read_uvh5(\n    files[0],\n    ant_str='cross',\n    times=times[0])\nfor file, time in zip(files[1:], times[1:]):\n    uvdi = UVData()\n    uvdi.read_uvh5(\n        file, \n        ant_str='cross',\n        times=time)\n    uvd += uvdi")


# In[6]:


get_ipython().run_cell_magic(u'time', u'', u'# Apply flags\nuvd.data_array *= np.logical_not(uvd.flag_array)')


# In[7]:


get_ipython().run_cell_magic(u'time', u'', u"# Intialize a cosmology and a beam\ncosmo = hp.conversions.Cosmo_Conversions()\nbeamfile = os.path.join(DATA_PATH, 'NF_HERA_Beams.beamfits')\nuvb = hp.pspecbeam.PSpecBeamUV(beamfile, cosmo=cosmo)")


# In[8]:


get_ipython().run_cell_magic(u'time', u'', u'# Convert to cosmological units (mK)\nJy_to_mK = uvb.Jy_to_mK(np.unique(uvd.freq_array), pol="xx")\nuvd.data_array *= Jy_to_mK[None, None, :, None]')


# In[9]:


get_ipython().run_cell_magic(u'time', u'', u'# Shift data and load datasets\nuvd1 = uvd.select(times=np.unique(uvd.time_array)[:-1:2], inplace=False)\nuvd2 = uvd.select(times=np.unique(uvd.time_array)[1::2], inplace=False)\nds = hp.PSpecData(dsets=[uvd1, uvd2], wgts=[None, None], beam=uvb)')


# In[10]:


get_ipython().run_cell_magic(u'time', u'', u"# Set visibility units\nds.dsets[0].vis_units = 'mK'\nds.dsets[1].vis_units = 'mK'")


# In[11]:


get_ipython().run_cell_magic(u'time', u'', u'# Phase data (What does this do?)\nds.rephase_to_dset(0)')


# In[121]:


get_ipython().run_cell_magic(u'time', u'', u'# Categorize baselines into physical separation length\nBIN_WIDTH = 0.3\nNORM_BINS = np.arange(0.0, 10000.0, BIN_WIDTH)\n\nantpos = {ant: pos for ant, pos in zip(uvd.get_ENU_antpos()[1], uvd.get_ENU_antpos()[0])}\n\nantpairs = uvd.get_antpairs()\n\nblpairs, blp_reds = [], {}\nbaselines, bls_reds = [], {}\nnorms = []\nfor antpair in antpairs:\n    ant0, ant1 = antpair\n    if (ant0 in xants) or (ant1 in xants) or (ant0 >= ant1):\n        continue\n    baselines.append(antpair)\n    blpair = (antpair, antpair)\n    blpairs.append(blpair)\n    norm = np.linalg.norm(antpos[ant0] - antpos[ant1])\n    norm = np.round(np.digitize(norm, NORM_BINS) * BIN_WIDTH, 1)\n    norms.append(norm)\n\n    if norm in bls_reds:\n        bls_reds[norm].append(antpair)\n        blp_reds[norm].append(blpair)\n    else:\n        bls_reds[norm] = [antpair]\n        blp_reds[norm] = [blpair]\nnorms = sorted(np.unique(norms))')


# In[124]:


blpairs


# In[13]:


get_ipython().run_cell_magic(u'time', u'', u'# Make UVPspec object\nuvp = ds.pspec(\n    baselines,\n    baselines,\n    (0, 1),\n    pols=("xx", "xx"),\n    spw_ranges=[(FREQrng[0], FREQrng[-1])],\n    taper="blackman-harris",\n    verbose=False)')


# In[14]:


get_ipython().run_cell_magic(u'time', u'', u'# Make UVPspec objects for each redundant baseline group\nuvps = []\nfor norm in norms:\n    uvpi = uvp.select(bls=bls_reds[norm], inplace=False)\n    uvps.append(uvpi)')


# In[115]:


# UVP = copy.deepcopy(uvp)
# BLP = copy.deepcopy(blp_reds.values())
# UVP = UVP.average_spectra(blpair_groups=BLP, inplace=False)

# BLPg = [bl[0] for bl in BLP]

# f, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(12, 12))
# hp.plot.delay_spectrum(
#     UVP,
#     BLPg,
#     0,
#     'xx',
#     average_blpairs=False,
#     average_times=True,
#     fold=True,
#     plot_noise=True,
#     delay=False,
#     deltasq=False,
#     legend=True,
#     ax=ax,
#     component='abs',
#     lines=True,
#     markers=False,
#     error=None)
# ax.legend(loc='upper right', ncol=3)


# In[110]:


# Average each UVPspec object in time and baseline bin and fold into wedge
for uvp, norm in zip(uvps, norms):
    blpairs = [[(bl, bl) for bl in bls_reds[norm]]]
    uvp.average_spectra(blpair_groups=blpairs, time_avg=True)
    uvp.fold_spectra()
    uvp.data_array[0] = uvp.data_array[0].reshape(
        (len(uvp.freq_array)))[np.nonzero(uvp.data_array[0].reshape((len(uvp.freq_array))))]


# In[111]:


# Load data from UVPspec objects into an array
wedge = np.array([uvp.data_array[0] for uvp in uvps])


# In[112]:


"""Plotting"""
def get_cmap(n, name='jet'):
    return plt.cm.get_cmap(name, n)
cmap = get_cmap(len(uvps))

# Make a copy of a uvp object to get freq_array and kparas
UVP = copy.deepcopy(uvps[0])

# Find bandwidth and central frequency in MHz for naming
BAND_START = (UVP.freq_array[0] * u.Hz).to(u.MHz)
BAND_STOP = (UVP.freq_array[-1] * u.Hz).to(u.MHz)
BANDWIDTH = (BAND_STOP - BAND_START)
CENTRAL_FREQ = ((BANDWIDTH / UVP.Nfreqs) + BAND_START)

# Generate x-values to plot against
kparas = (UVP.get_kparas(0)/u.Mpc).insert(0, 0)

Tsys = 360
plt.figure(figsize=(10, 10))
for i, (uvp, norm, pspec) in enumerate(zip(uvps, norms, wedge)):
    plt.plot(
        kparas,
        np.log10(np.abs(pspec)),
        c=cmap(i),
        ls='-',
        lw=1,
        label='{norm}m ({bls} bls)'.format(norm=norm, bls=len(bls_reds[norm])))

    noise = uvp.generate_noise_spectra(0, 'xx', Tsys)
    noise = noise[noise.keys()[0]]
    noise = np.insert(noise, 0, noise[0, 0], axis=1)
    noise = noise.reshape(len(kparas))
    noise = np.log10(noise)
    plt.plot(
        kparas,
        noise,
        c=cmap(i),
        ls='--',
        lw=1,
        label='{norm}m {Tsys}K'.format(norm=norm, Tsys=Tsys))
    
    
plt.legend(loc='upper right', ncol=3)

# x-axis
plt.xlim((0, UVP.get_kparas(0)[-1]))
plt.xlabel(r"$k_{\parallel}\ [\rm\ Mpc^{-1}\ h]$", size=20)

# y-axis
plt.ylim((8, 15))
plt.ylabel(r"$P(k)\ \rm [\log_{10}({mK^2\ Mpc^3\ h^{-3}})]$", size=20)

# Titles
plt.title("XX")
plt.suptitle("Bandwidth: {BW}\nCentral Frequency: {CF}\nExcluded Antennae: {XANTS}".format(
    BW=np.round(BANDWIDTH, 2),
    CF=np.round(CENTRAL_FREQ, 1),
    XANTS=', '.join([str(xant) for xant in xants])))

# Save and show the plot with a grid
plt.grid()
plt.savefig("zen.{JD}.{JD0}_{JDf}.xx.HH{DFext}.pdf".format(
    BW=np.round(BANDWIDTH.value, 2),
    CF=np.round(CENTRAL_FREQ.value, 1),
    JD=files[0].split(".")[1],
    JD0=files[0].split(".")[2],
    JDf=files[-1].split(".")[2],
    DFext=os.path.splitext(files[0])[1]))
plt.show()


# In[ ]:


# plt.figure(figsize=(10, 10))
# plt.imshow(np.log10(np.abs(wedge)), interpolation="nearest", aspect="auto")

# plt.tick_params(axis='both', direction='inout')

# plt.xticks([])
# plt.xlabel(str(norms[i]) + " m", rotation=45, ha="center")

# horizon = ((norms[i]*u.m / c.c).to(u.ns)).value
# plt.axhline(y=horizon, color="w", ls=":")
# plt.axhline(y=-horizon, color="w", ls=":")

# plt.ylim((uvp.get_dlys(0)[0]*1e9 / 2., uvp.get_dlys(0)[-1]*1e9 / 2.))
    
# plt.text(0.07, 0.5, r"$\tau$ [ns]", ha="center", rotation="vertical", size=20)
# plt.text(0.5, 0.04, "Redundant Baseline Group", ha="center", size=20)
# plt.subplots_adjust(wspace=0, hspace=0)
# cbar_ax = fig.add_axes([0.9125, 0.25, 0.025, 0.5])
# cbar = fig.colorbar(im, cax=cbar_ax)
# cbar.set_label(r"$P(k)\ \rm [mK^2\ h^{-3}\ Mpc^3]$", fontsize=20, ha='center')
# plt.show()
# plt.clf()


# In[ ]:


# Plot each baselines divided by or subtracted by its noise estimate in 2d wedge format

