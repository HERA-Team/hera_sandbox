import aipy as a
import sys
import capo
from glob import glob
import pylab as plt
import numpy as np
import time

def flag_avg(data,flags):
    avg = (data*flags).sum(axis=0)/flags.sum(axis=0)
    return avg

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

folio = False
if folio:
    datadir = '/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/'
    outdir = '/data4/paper/jaguirre/GlobalSignal/'
else:
    datadir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/Omnicaled/'
    outdir = '/Users/jaguirre/Documents/PAPER/GlobalSignal/BaselinePulls/'
    
files = sorted(glob(datadir+'zen.2456680.[2,3,4,5]*.uvcRREO'))

ants = np.arange(0,128)
nant = len(ants)
nbl = nant*(nant+1)/2.
ibl = 1
polstr = 'xx'

anttuples = []
for i in ants:
    for j in np.arange(i,nant):
        anttuples.append((i,j))
antstrings = tuplist2strlist(anttuples)

#%%
avgspecs={}
nbl_per_chunk = 64
nchunks = int(nbl/nbl_per_chunk)
for i in np.arange(nchunks):
    chunk = np.arange(i*nbl_per_chunk,(i+1)*nbl_per_chunk)
    antstr = ''
    for c in chunk:
        antstr+=antstrings[c]+','
    antstr = antstr[0:-1]   
    t0 = stime(message='Getting '+antstr+'; chunk '+str(i+1)+' of '+str(nchunks))
    tinfo,wfall,flags = capo.arp.get_dict_of_uv_data(files,antstr=antstr,polstr=polstr)
    etime(t0)
    t0 = stime(message='Writing out')
    for c in chunk:
        anttup = anttuples[c]
        w = wfall[anttup][polstr]
        ifl = np.ones(flags[anttup][polstr].shape)
        ifl[flags[anttup][polstr]] = 0
        avgspec = flag_avg(w,ifl)
        avgspecs[anttup] = avgspec
        fs = antstrings[c].split('_')
        filestr = fs[0].zfill(3)+'_'+fs[1].zfill(3)
        np.savez(outdir+'bl_'+filestr+'.npz',wfall=w,iflags=ifl,avgspec=avgspec)
        np.savez(outdir+'avgspecs.npz',avgspecs=avgspecs)
    etime(t0)
    #w = {}
    #iflags = {}
    #anttuple = anttuples[chunk]
    #for anttup in anttuple:
    #w[anttup] = wfall[anttup]['xx']
    #ifl = np.ones(flags[anttup]['xx'].shape)
    #ifl[flags[anttup]['xx']] = 0
    #iflags[anttup] = ifl
    
