#! /usr/bin/env python
import os
import argparse
import pyuvdata
import numpy as np
import pyuvdata.utils as uvutils
import matplotlib.pyplot as plt

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

plt.scatter(antpos[:, 0], antpos[:, 1], c=amps)
plt.colorbar()
for ant in range(uv.Nants_telescope):
        plt.annotate(str(uv.antenna_numbers[ant]), xy=antpos[ant, 0:2], textcoords='data')
plt.show()
