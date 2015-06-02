#! /usr/bin/env python

###
# NAME:    check_baselines.py
# AUTHOR:  Carina Cheng
# PURPOSE: Extracts a certain type of baseline from PSA128 slice data (npz files) and compares LST models between them
###

import numpy
import matplotlib.pyplot as plt
import aipy
import os, sys
import optparse

o = optparse.OptionParser()
o.set_usage('check_baselines.py [options] *npz')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

firstfile = numpy.load(args[0])
keys = firstfile.keys()
bls = [] #baselines by aipy convention
bls_ants = [] #baselines by antenna numbers

for i in range(len(keys)):
    if keys[i][0] != 't':
        bls.append(keys[i])
        bls_ants.append(aipy.miriad.bl2ij(keys[i]))

### All baselines in the npz slices are 30-m E/W baselines ###
### There are 98 of them ###

mybl = str(aipy.miriad.ij2bl(64,49))
all_data = []
all_lsts = []

for filename in args:
    jd_int = int(filename.split('.')[1])
    jd_dec = float(filename.split('.')[1]+'.'+filename.split('.')[2])
    npz = numpy.load(filename)
    try:
        data = npz[mybl]
        lsts = npz['t'+mybl]*12/numpy.pi #LST hours  
        all_data.append(data)
        all_lsts.append(lsts)
    except:
        continue

all_data = numpy.array(all_data)
all_lsts = numpy.array(all_lsts)
all_data = numpy.concatenate(all_data)
all_lsts = numpy.concatenate(all_lsts)

#binning LST
all_lsts,all_data = zip(*sorted(zip(all_lsts,all_data)))
num_bins = 100
d_lst = float((all_lsts[-1]-all_lsts[0])/num_bins)
lst_grid = numpy.arange(all_lsts[0],all_lsts[-1]+(d_lst/2),d_lst)
model = numpy.zeros(num_bins)
for i in range(num_bins):
    bin = []
    for j in range(len(all_lsts)):
        if all_lsts[j] >= lst_grid[i] and all_lsts[j] < lst_grid[i+1]:
            bin.append(all_data[j])
    med_val = numpy.median(bin)
    model[i] = med_val.real

#plotting LST model
plt.plot(lst_grid[:-1],model,'k-')
plt.show()

#TO DO:
   # loop through baselines
   # have command-line option to designate 'good' baseline
   # plot all baseline models on same plot?
   # quantitative way to figure out which baselines deviate from the 'good' baseline?


