# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 14:21:37 2016

@author: jaguirre
"""

import aipy
import sys
import capo
from glob import glob
import pylab as plt
import numpy as np
import time
import h5py
#%%
def stime(message=''):
    print message
    t0 = time.time()
    return t0

def etime(t0):
    dt = time.time() - t0
    print 'Took',dt,'seconds'
    return dt

caldir='/Users/jaguirre/PyModules/capo/jea/cals/'
sys.path.append(caldir)
calfile = 'psa6622_v003'
bldir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/BaselinePulls/'
hdf5dir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/'
#%%

ants = np.arange(0,128)
nant = len(ants)
nbl = nant*(nant+1)/2.
polstr = 'xx'
nfreqs = 203

freqs = np.linspace(100.,200.,num=nfreqs)

aa = aipy.cal.get_aa(calfile,np.array([.15]))
info = capo.omni.aa_to_info(aa)
reds = info.get_reds()
#%%
#xnom = info.antloc[:,1]*15.
#ynom = info.antloc[:,0]*4.
#plt.figure(1)
#plt.clf()
#plt.plot(xnom,ynom,'o')
#plt.xlim([-15,240])
#plt.ylim([-4,30])
#%%
ibl = 1
f = h5py.File(hdf5dir+'zen.2456680.xx.uvcRREO.hdf5','w')
wfalls = f.create_group("waterfalls")
flags = f.create_group("iflags")

avgspecs = np.zeros([nant,nant,203],dtype='complex128')
for i in ants:
    for j in np.arange(i,nant):
        print 'Baseline',ibl,'of',nbl
        blstr = (str(i).zfill(3)+'_'+str(j).zfill(3))
        tup = (i,j)
        data = np.load(bldir+'bl_'+blstr+'.npz')
        wfalls.create_dataset(blstr,data=data['wfall'])
        flags.create_dataset(blstr,data=data['iflags'])        
        avgspecs[i,j,:]=data['avgspec']
        avgspecs[j,i,:]=data['avgspec']
        ibl+=1
        
f.close()
# This got fouled up
#avgspecs=np.load(bldir+'avgspecs.npz')
#%%
# Find (10,64)
myblgrp = -1
for i,group in enumerate(reds):
    for bl in group:
        if (bl == (10,64) or bl == (64,10)):
            myblgrp = i
#%%
badants = [18,5,61,72,76,77,78,79,16,26,34,38,46,50,84,99]
plt.figure(1)
plt.clf()
for bl in reds[myblgrp]:
    i = bl[0]
    j = bl[1]
    if i in badants or j in badants:
        continue
    plt.plot(freqs,np.abs(avgspecs[i,j,:]))
    


    