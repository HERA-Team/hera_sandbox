#! /usr/bin/env python
import os
import argparse
import pyuvdata
import numpy as np
import pyuvdata.utils as uvutils
import matplotlib.pyplot as plt
from astropy import time
from hera_mc import corr_handling
import re

parser = argparse.ArgumentParser(description='Plot auto locations and magnitudes')
parser.add_argument('file', metavar='file', type=str, nargs='+',
                    help='File to be processed.')

args = parser.parse_args()

uv = pyuvdata.UVData()
uv.read_miriad(args.file)
antpos = uv.antenna_positions + uv.telescope_location
antpos = uvutils.ENU_from_ECEF(antpos.T, *uv.telescope_location_lat_lon_alt).T

amps = np.zeros(uv.Nants_telescope)
for ant in range(uv.Nants_telescope):
    d = uv.get_data((uv.antenna_numbers[ant], uv.antenna_numbers[ant]))
    amps[ant] = np.mean(np.abs(d))

at_time = time.Time(uv.extra_keywords['obsid'], format='gps')
h = corr_handling.Handling()
pol = uvutils.polnum2str(uv.polarization_array[0])[0]
if pol == 'X':
    pol = 'e'
else:
    pol = 'n'

f = plt.figure(figsize=(10, 8))
plt.scatter(antpos[:, 0], antpos[:, 1], c=amps)
plt.clim([0, amps.max()])
plt.colorbar()
receiverators = []
pams = []
texts = []
for ant in range(uv.Nants_telescope):
    pam = h.get_pam_info(uv.antenna_names[ant], at_time)
    text = (str(uv.antenna_numbers[ant]) + pol + '\n' + pam[pol][0] + '\n' +
            pam[pol][1])
    texts.append(text)
    result = re.findall(r'RI(\d+)', pam[pol][0])[0]
    receiverators.append(int(result))
    pams.append(pam[pol][1])
    plt.annotate(text, xy=antpos[ant, 0:2] + [1, 0], textcoords='data')
xr = antpos[:, 0].max() - antpos[:, 0].min()
yr = antpos[:, 1].max() - antpos[:, 1].min()
plt.xlim([antpos[:, 0].min() - 0.05 * xr, antpos[:, 0].max() + 0.2 * xr])
plt.ylim([antpos[:, 1].min() - 0.05 * yr, antpos[:, 1].max() + 0.1 * yr])
plt.title(str(at_time.datetime))

f, axarr = plt.subplots(1, np.max(receiverators), sharex=True, sharey=True, figsize=(15, 7))
for rxr in np.unique(receiverators):
    ind = np.where(receiverators == rxr)[0]
    for i, ant in enumerate(ind):
        ax = axarr[rxr - 1].scatter(0, i, c=amps[ant], vmin=0, vmax=amps.max())
        axarr[rxr - 1].annotate(str(uv.antenna_numbers[ant]) + ',' + pams[ant], xy=[0.01, i])
    axarr[rxr - 1].set_yticks([])
    axarr[rxr - 1].set_xticks([])
for rxr in range(np.max(receiverators)):
    axarr[rxr].set_title('Rxr ' + str(rxr + 1))
plt.xlim([-.01, .1])
plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(wspace=0)
cbar_ax = f.add_axes([.14, .05, .72, .05])
f.colorbar(ax, cax=cbar_ax, orientation='horizontal')
f.suptitle(str(at_time.datetime))
plt.show()
