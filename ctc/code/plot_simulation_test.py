#!/usr/bin/env python

import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt
import optparse
import os, sys

#plotting

plt.figure(figsize = (12,8))
#plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.3)


file1 = aipy.miriad.UV('/Users/carinacheng/capo/ctc/tables/vis_simulation_0630_v1.uv')
file2 = aipy.miriad.UV('/Users/carinacheng/capo/ctc/tables/vis_simulation_0702_v2.uv')
file3 = aipy.miriad.UV('/Users/carinacheng/capo/ctc/tables/vis_simulation_0702_v1.uv')

data1 = []
times1 = []
for p,d,f in file1.all(raw=True):
    data1.append(d[0])
    times1.append(p[1])

data2 = []
times2 = []
for p,d,f in file2.all(raw=True):
    data2.append(d[30:170]) 
    times2.append(p[1])

data3 = []
times3 = []
for p,d,f in file3.all(raw=True):
    data3.append(d[30:170]) 
    times3.append(p[1])

#plt.subplot(2,1,1)
p1 ,= plt.plot(times1,numpy.real(data1),'r-',label='test')
plt.plot(times2,numpy.real(data2),'k-')
plt.plot(times3,numpy.real(data3),'b-')
#plt.ylim(-1e-6,1e-6)
plt.xlim(2454500,2454500.5)
plt.xlabel("Time (Julian Date)")
plt.ylabel("Visibilities")
#plt.legend([p1,p2,p3],['f=0.15GHz','all f, no interp','all f, interp'], prop={'size':11},loc=2)
plt.show()

