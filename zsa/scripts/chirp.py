#! /usr/bin/env python

import numpy
import aipy, capo
import matplotlib.pyplot as plt
import glob
import sys, os

###
# Investigates the chirp found in S1E2 data
# Plots the gains of all antennas for a time integration when the chirp is present
###


#Read files and load data
#files = numpy.sort(glob.glob('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/zen.2456684.[4-5]*xx.npz')) #XXX for these files, the chirp occurs prominently at times 0-150 
files = numpy.sort(glob.glob('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_test_ctc/zen.2456722.4*.xx.npz'))
#files = numpy.sort(glob.glob('/home/zakiali/omni_outputs/zen.2457458.2*npz'))

ants = capo.omni.from_npz(str(files[0]))[1]['x'].keys()
gains = {}
vis = {}
xtalk = {}
for f in files:
    print 'Reading', f
    m,g,v,x = capo.omni.from_npz(str(f))
    for a in g['x'].keys():
        try:
            gains[a] = numpy.vstack((gains[a],g['x'][a]))
        except: gains[a] = g['x'][a]
    for bl in v['xx'].keys():
        try:
            vis[bl] = numpy.vstack((vis[bl],v['xx'][bl]))
            xtalk[bl] = numpy.vstack((xtalk[bl],x['xx'][bl]))
        except: 
            vis[bl] = v['xx'][bl]
            xtalk[bl] = x['xx'][bl]

ant1 = 1
ant2 = 4
plt.clf()
plt.subplot(1,4,1)
plt.imshow(numpy.real(gains[ant1]),aspect='auto')
#plt.imshow(numpy.angle(gains[ant1]),aspect='auto')
plt.title('real(gain) for antenna '+str(ant1))
plt.subplot(1,4,2)
#plt.imshow(numpy.angle(vis[(ant1,ant2)]),aspect='auto')
capo.arp.waterfall(vis[(ant1,ant2)],mode='log',drng=4,mx=0)
plt.title('mdlvis for ('+str(ant1)+','+str(ant2)+')')
plt.subplot(1,4,3)
#plt.imshow(numpy.angle(gains[ant1]*numpy.conj(gains[ant2])*vis[(ant1,ant2)]),aspect='auto')
capo.arp.waterfall(gains[ant1]*numpy.conj(gains[ant2])*vis[(ant1,ant2)],mode='log',drng=4,mx=0)
plt.title('g'+str(ant1)+'g'+str(ant2)+'*mdlvis')
plt.subplot(1,4,4)
capo.arp.waterfall(xtalk[(ant1,ant2)],mode='log',drng=4,mx=0)
plt.title('xtalk for ('+str(ant1)+','+str(ant2)+')')
plt.show()




"""
#Subtract off background
for a in gains.keys():
    median = numpy.median(gains[a],axis=0)
    gains[a] -= median #subtract off median in time (background)
    gains[a] = numpy.real(gains[a])
    mk = numpy.ma.masked_where(gains[a] == 1,gains[a]).mask #flags
    gains[a] = numpy.ma.masked_array(gains[a],mask=mk) #masked array

#Plot background-subtracted gains for all antennas
subplotnum = 1
plotnum = 1
plt.figure(plotnum,figsize=(10,10))
for a in gains.keys():
    if subplotnum == 26:
        #break #only generate one page of plots (faster for testing) 
        plotnum += 1
        plt.figure(plotnum,figsize=(10,10))
        subplotnum = 1
    plt.subplot(5,5,subplotnum)
    plt.imshow(gains[a],aspect='auto',interpolation='nearest')
    plt.title(a,fontsize=10)
    plt.tick_params(axis='both',which='major',labelsize=6)
    plt.tight_layout()
    subplotnum += 1
plt.show()
plt.clf()
"""
"""
#Plot slice in time for all antennas

plt.imshow(gains[0],aspect='auto')
plt.title('real(gain)-median in time for antenna 0')
plt.show()
for a in gains.keys():
    plt.plot(gains[a][95],label=a)
plt.legend(prop={'size':8})
plt.show()
"""




