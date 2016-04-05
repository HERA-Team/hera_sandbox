#! /usr/bin/env python

import numpy
import aipy
import capo
import glob
import matplotlib.pyplot as plt

#Tests a simple xtalk subtraction (subtracts the average of a few baselines)
#Used to see why Omnical doesn't seem to be working


files = numpy.sort(glob.glob('/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v1_xtalk/zen.2456942.*xx.uvcRREOA'))

print 'Reading', str(len(files)),'files...'

t,d,f = capo.arp.get_dict_of_uv_data(files,antstr='49_64,1_4,105_106,11_36,1_58,66_67,9_58,14_54,102_103',polstr='xx')
#d_6667 = d[aipy.miriad.ij2bl(66,67)]['xx']
#d[aipy.miriad.ij2bl(1,4)]['xx'] -= d_6667 #hacking script to look at diff between two baselines

alldata = []

plt.figure(1)
for i,bl in enumerate(d.keys()):
    plt.subplot(3,3,i+1)
    capo.arp.waterfall(d[bl]['xx'],mx=0,drng=4)
    alldata.append(d[bl]['xx'])
    plt.colorbar()
    plt.title(aipy.miriad.bl2ij(bl))
plt.suptitle('Omnical_v1 Results')

print 'Performing Manual xtalk Removal...'

alldata = numpy.array(alldata)
avg = numpy.mean(alldata,axis=0) #avg for specified baselines

plt.figure(2)
for i,bl in enumerate(d.keys()):
    plt.subplot(3,3,i+1)
    res = numpy.mean(d[bl]['xx']-avg,axis=0) #averaged over all of time (not per file unless only one file is specified)
    res = numpy.resize(res,(d[bl]['xx'].shape[0],res.shape[0]))
    capo.arp.waterfall(d[bl]['xx']-res,mx=0,drng=4)
    plt.colorbar()
    plt.title(aipy.miriad.bl2ij(bl))
plt.suptitle('Omnical_v1 - Averaged Residual (in time)')

plt.show()
