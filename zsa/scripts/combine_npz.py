#! /usr/bin/env python
import numpy as np
import sys, aipy as a

files = sys.argv[1:]

badants = [2,10,14,15,16,22,26,27,28,31,33,34,38,42,43,46,47,50,53,58,64,72,74,84,91,97,105,107] # epoch1 bad antennas
#goodjds = np.loadtxt('/home/zakiali/hackweek/good_days_epoch1_xx.txt')
goodjds = np.loadtxt('/home/zakiali/hackweek/good_days_epoch2.txt')
CHANS = 40
print goodjds
samples = {}
for f in files:
    cmpstr = f.split('.')[1]
    if not float(cmpstr) in goodjds : continue
    if f.split('.')[2].startswith('17'): print 'reading ',f 
    d = np.load(f)
    for ant in d.keys():
        if ant.isdigit():
            if not ant in samples.keys(): samples[ant] = d[ant][:,:CHANS]
            else: samples[ant] = np.concatenate([samples[ant],d[ant][:,:CHANS]])

for ant in samples.keys():
    hreadable = a.miriad.bl2ij(ant)[0]
    times = samples[ant].shape[0]
    d,m = divmod(times,CHANS)
    gooddir='good2/'
    baddir='bad2/'
    if a.miriad.bl2ij(ant)[0] in badants:
        np.savez(baddir+'%d_data.npz'%hreadable, samples[ant][:d*CHANS,:].reshape(-1,CHANS,CHANS))
    else: np.savez(gooddir+'%d_data.npz'%hreadable, samples[ant][:d*CHANS,:].reshape(-1,CHANS,CHANS))

