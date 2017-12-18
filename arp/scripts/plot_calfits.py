#! /usr/bin/env python
import numpy as np, pylab as plt
import pyuvdata, hera_cal
import uvtools
import sys

omni_files = sys.argv[1:]
print omni_files

oc = pyuvdata.UVCal()
oc.read_calfits(omni_files)
gains = oc.gain_array.squeeze()
print gains.shape
for i in [0,1,5]:
    ai = oc.ant_array[i]
    g = gains[i,...,0].T
    plt.figure()
    uvtools.plot.waterfall(g, mode='lin', mx=2, drng=2)
    plt.title(ai)
    plt.figure()
    uvtools.plot.waterfall(g, mode='phs')
    plt.title(ai)
    plt.show()
