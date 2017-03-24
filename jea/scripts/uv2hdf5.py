#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:16:39 2016

@author: jaguirre
"""

import aipy as a
import sys
import capo
from glob import glob
import pylab as plt
import numpy as np
import time
import h5py

def tuplist2str(tuplist):
    string = ''
    for tup in tuplist:
        string += str(tup[0])+'_'+str(tup[1])+','
    string = string[0:-1]
    return string

def tuplist2strlist(tuplist):
    stringlist= []
    for tup in tuplist:
        stringlist.append(str(tup[0])+'_'+str(tup[1]))
    return stringlist

def stime(message=''):
    print message
    t0 = time.time()
    return t0

def etime(t0):
    dt = time.time() - t0
    print 'Took',dt,'seconds'
    return dt

# Define the data file    
datadir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/RFIdf6fchan_HH/'
outdir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/'
files = sorted(glob(datadir+'*.uvcR'))

#%%
# Define the antennas
ants = np.array([9, 10, 20, 22, 31, 43, 53, 64, 65, 72, 80, 81, 88, 89, 96, 97, 104, 105, 112])#np.arange(0,128)
nant = len(ants)
nbl = nant*(nant+1)/2.
ibl = 1
polstr = 'xx'

anttuples = []
for ni,i in enumerate(ants):
    for j in ants[ni:]: #np.arange(i,nant):
        anttuples.append((i,j))
antstrings = tuplist2strlist(anttuples)

#%% Antenna information
#nfreqs = 203
#freqs = np.linspace(100.,200.,num=nfreqs)
#
#aa = aipy.cal.get_aa(calfile,np.array([.15]))
#info = capo.omni.aa_to_info(aa)
#reds = info.get_reds()

#%%
#avgspecs={}
nbl_per_chunk = 64 # This is tuned for speed & memory on a given processor 
nchunks = int(nbl/nbl_per_chunk)

# Open the hdf file and create the groups
f = h5py.File(outdir+'HERA.hdf5','w')
wfalls = f.create_group("waterfalls")
flags = f.create_group("iflags")

for i in np.arange(nchunks):
    chunk = np.arange(i*nbl_per_chunk,(i+1)*nbl_per_chunk)
    antstr = ''
    for c in chunk:
        antstr+=antstrings[c]+','
    antstr = antstr[0:-1]   
    t0 = stime(message='Getting '+antstr+'; chunk '+str(i+1)+' of '+str(nchunks))
    tinfo,wfall,flag = capo.arp.get_dict_of_uv_data(files,antstr=antstr,polstr=polstr)
    etime(t0)
    t0 = stime(message='Writing out')
    for c in chunk:
        anttup = anttuples[c]
        w = wfall[anttup][polstr]
        fl = flag[anttup][polstr]
        print 'anttup', anttup
        ant_i = anttup[0]
        ant_j = anttup[1]
        blstr = (str(ant_i).zfill(3)+'_'+str(ant_j).zfill(3))
        wfalls.create_dataset(blstr,data=w)
        flags.create_dataset(blstr,data=fl)  
        ibl += 1
    etime(t0)
    
f.close()
