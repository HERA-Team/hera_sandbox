#! /usr/bin/env python

import aipy
import numpy
import capo
import glob
import matplotlib.pyplot as plt
import os, sys

###
# Compares phases of Omnicaled data from version 1 vs. version 2
###

#Epoch 3 Comparison Day: JD 2456962
v2path = numpy.sort(glob.glob('/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2/zen.2456942.[3-4]*.xx.uvcRREO'))
v1path = numpy.sort(glob.glob('/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/good/zen.2456942.[3-4]*.xx.uvcRREOO'))
print 'Reading', v1path
print 'Reading', v2path

#get data
t1,d1,f1 = capo.arp.get_dict_of_uv_data(v1path,antstr='64_49',polstr='xx',return_lsts=True)
t2,d2,f2 = capo.arp.get_dict_of_uv_data(v2path,antstr='64_49',polstr='xx',return_lsts=True)
d1 = d1[aipy.miriad.ij2bl(64,49)]['xx']
d2 = d2[aipy.miriad.ij2bl(64,49)]['xx']

print d1.shape,d2.shape

#plot data
plt.subplot(221)
capo.arp.waterfall(d1,drng=3,mode='phs',mx=0)
plt.title('Omni_v1 Phase')
plt.subplot(222)
capo.arp.waterfall(d2,drng=3,mode='phs',mx=0)
plt.title('Omni_v2 Phase')
plt.subplot(223)
capo.arp.waterfall(d1,drng=3,mode='log',mx=0)
plt.title('Omni_v1 Waterfall')
plt.subplot(224)
capo.arp.waterfall(d2,drng=3,mode='log',mx=0)
plt.title('Omni_v2 Waterfall')
plt.tight_layout()
plt.show()

#get phases
d1 = numpy.angle(d1)
d2 = numpy.angle(d2)

#plot phase differences
plt.imshow(d1-d2)
plt.colorbar()
plt.title('Omni_v1 Phases - Omni_v2 Phases')
plt.show()
