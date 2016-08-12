#! /usr/bin/env python

import numpy
import aipy, capo
import matplotlib.pyplot as plt
import glob
import optparse
import sys, os

###
# Investigates Omnical .npz outputs
# Plots the gain of antenna, model visibility of a baseline, gains applied to the model visibility, and xtalk
###

o = optparse.OptionParser()
o.set_usage('plot_npz.py *npz')
aipy.scripting.add_standard_options(o,ant=True,pol=True)
o.set_description(__doc__)

opts,args = o.parse_args(sys.argv[1:])

#Read files and load data
#files = numpy.sort(glob.glob('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/zen.2456684.[4-5]*xx.npz')) #XXX for these files, the chirp occurs prominently at times 0-150 
#files = numpy.sort(glob.glob('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/zen.2456722.4*.xx.npz')

#Get data
m,g,v,x = capo.omni.from_npz(args,verbose=True)
ant1,ant2 = opts.ant.split('_')
ant1 = int(ant1)
ant2 = int(ant2)
pol = opts.pol

#Plot
plt.clf()
plt.subplot(2,4,1)
plt.imshow(numpy.real(g[pol[0]][ant1]),aspect='auto')
plt.title('real(gain) for antenna '+str(ant1))
plt.subplot(2,4,2)
capo.arp.waterfall(v[pol][(ant1,ant2)],mode='log',drng=4,mx=0)
plt.title('mdlvis for ('+str(ant1)+','+str(ant2)+')')
plt.subplot(2,4,3)
capo.arp.waterfall(g[pol[0]][ant1]*numpy.conj(g[pol[0]][ant2])*v[pol][(ant1,ant2)],mode='log',drng=4,mx=0)
plt.title('g'+str(ant1)+'g'+str(ant2)+'*mdlvis')
plt.subplot(2,4,4)
capo.arp.waterfall(x[pol][(ant1,ant2)],mode='log',drng=4,mx=0)
plt.title('xtalk for ('+str(ant1)+','+str(ant2)+')')

plt.subplot(2,4,5)
plt.imshow(numpy.angle(g[pol[0]][ant1]),aspect='auto')
plt.title('phase of gains for ant %d'%ant1)
plt.subplot(2,4,6)
plt.imshow(numpy.angle(v[pol][(ant1,ant2)]),aspect='auto')
plt.title('phase of mdlvis')
plt.subplot(2,4,7)
plt.imshow(numpy.angle(g[pol[0]][ant1]*numpy.conj(g[pol[0]][ant2])*v[pol][(ant1,ant2)]),aspect='auto')
plt.title('phase of g'+str(ant1)+'g'+str(ant2)+'*mdlvis')
plt.subplot(2,4,8)
plt.imshow(numpy.angle(x[pol][(ant1,ant2)]),aspect='auto')
plt.title('phase of xtalk')

plt.show()




