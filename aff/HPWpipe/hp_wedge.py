
# coding: utf-8

# In[1]:


import os
import copy
import glob

import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt

import hera_pspec as hp
from hera_pspec.data import DATA_PATH
from pyuvdata import UVData
from hera_cal import redcal


# In[2]:


# Locate data
dfiles = sorted(glob.glob("/lustre/aoc/projects/hera/H1C_IDR2/IDR2_1/2458111/zen.2458111.30???.xx.HH.uvOC"))


# In[3]:


get_ipython().run_cell_magic('time', '', '# Initialize UVData object to read and contain MIRIAD information\nuvd = UVData()\nuvd.read_miriad(dfiles, ant_str=\'cross\')\n\n# Apply flags\nuvd.data_array *= np.logical_not(uvd.flag_array)\n\n# Intialize a cosmology and a beam\ncosmo = hp.conversions.Cosmo_Conversions()\nbeamfile = os.path.join(DATA_PATH, \'NF_HERA_Beams.beamfits\')\nuvb = hp.pspecbeam.PSpecBeamUV(beamfile, cosmo=cosmo)\n\n# Convert to cosmological units (mK)\nJy_to_mK = uvb.Jy_to_mK(np.unique(uvd.freq_array), pol="xx")\nuvd.data_array *= Jy_to_mK[None, None, :, None]\n\n# Shift data and load datasets\nuvd1 = uvd.select(times=np.unique(uvd.time_array)[:-1:2], inplace=False)\nuvd2 = uvd.select(times=np.unique(uvd.time_array)[1::2], inplace=False)\nds = hp.PSpecData(dsets=[uvd1, uvd2], wgts=[None, None], beam=uvb)\n\nds.dsets[0].vis_units = \'mK\'\nds.dsets[1].vis_units = \'mK\'')


# In[4]:


get_ipython().run_cell_magic('time', '', '# Phase data (What does this do?)\nds.rephase_to_dset(0)')


# In[5]:


get_ipython().run_cell_magic('time', '', '# Categorize baselines into physical separation length\nBIN_WIDTH = 0.3\nNORM_BINS = np.arange(0.0, 10000.0, BIN_WIDTH)\n\nantpos = {ant: pos for ant, pos in zip(uvd.get_ENU_antpos()[1], uvd.get_ENU_antpos()[0])}\n\nantpairs = uvd.get_antpairs()\nxants1, xants2 = hp.utils.calc_reds(uvd1, uvd2)[3:]\nxants = np.unique(xants1 + xants2)\n\nreds = {}\nfor antpair in antpairs:\n    ant0, ant1 = antpair\n    if (ant0 in xants) or (ant1 in xants) or (ant0 >= ant1):\n        continue\n    norm = np.linalg.norm(antpos[ant0] - antpos[ant1])\n    norm = np.round(np.digitize(norm, NORM_BINS) * BIN_WIDTH, 1)\n\n    if norm in reds:\n        reds[norm].append(antpair)\n    else:\n        reds[norm] = [antpair]\nnorms = sorted(reds.keys())')


# In[6]:


get_ipython().run_cell_magic('time', '', '# Initialize UVPspec objects for each baseline bin\nuvps = [] \nfor norm in norms:\n    uvp = ds.pspec(\n        reds[norm],\n        reds[norm],\n        (0, 1),\n        pols=("xx", "xx"),\n        spw_ranges=[(580,680)],\n        taper="blackman-harris",\n        verbose=False)\n    uvps.append(uvp)')


# In[51]:


np.savez(
    'metadata.npz',
    reds=reds,
    norms=norms,
    xants=xants,
    antpos=antpos)


# In[ ]:


get_ipython().run_cell_magic('time', '', '# Average each UVPspec object in time and baseline bin and fold into wedge\nfor uvp, norm in zip(uvps, norms):\n    blpairs = [[(bl, bl) for bl in reds[norm]]]\n    uvp.average_spectra(blpair_groups=blpairs, time_avg=True)\n    uvp.fold_spectra()\n    uvp.data_array[0] = uvp.data_array[0].reshape(\n        (len(uvp.freq_array)))[np.nonzero(uvp.data_array[0].reshape((len(uvp.freq_array))))]')


# In[ ]:


# Load data from UVPspec objects into an array
wedge = np.array([uvp.data_array[0] for uvp in uvps])


# In[ ]:


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

Tsys = 400
plt.figure(figsize=(10, 10))
for i, (uvp, norm, pspec) in enumerate(zip(uvps, norms, wedge)):
    plt.plot(
        kparas,
        np.log10(np.abs(pspec)),
        c=cmap(i),
        ls='-',
        lw=1,
        label='{norm}m ({ants} ants)'.format(norm=norm, ants=len(reds[norm])))

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
# plt.ylim((0, 20))
plt.ylabel(r"$P(k)\ \rm [\log_{10}({mK^2\ Mpc^3\ h^{-3}})]$", size=20)

# Titles
plt.title("pol: xx; Bandwidth: {BW}; Central Frequency: {CF}".format(
    BW=np.round(BANDWIDTH, 2),
    CF=np.round(CENTRAL_FREQ, 1)))
plt.suptitle(os.path.basename(dfiles[0]) + "\nto\n" + os.path.basename(dfiles[-1]))

# Save and show the plot with a grid
plt.grid()
plt.savefig("xx.{BW}_{CF}.{JD0}_{JDf}{DFext}.pdf".format(
    BW=np.round(BANDWIDTH.value, 2),
    CF=np.round(CENTRAL_FREQ.value, 1),
    JD0=dfiles[0].split(".")[2],
    JDf=dfiles[-1].split(".")[2],
    DFext=os.path.splitext(dfiles[0])[1]))
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

