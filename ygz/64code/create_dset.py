#! /usr/bin/env python

import numpy as np, aipy, capo, pylab as plt, sys, glob
import md5
from joblib import Parallel, delayed
import multiprocessing, joblib
import oqe
import os
from waterfall import waterfall
num_cores = multiprocessing.cpu_count()
def dB(sig): return 10*np.log10(np.abs(np.average(sig.real, axis=1)))

aa = aipy.cal.get_aa('psa6622_v001',np.array([.15]))
def find_sep(aa, bls, drow=None, dcol=None):
    layout = aa.ant_layout
    rv = []
    for i,j in bls:
        irow,icol = np.where(layout == i)
        jrow,jcol = np.where(layout == j)
        if not drow is None and abs(irow - jrow) != drow: continue
        if not dcol is None and abs(icol - jcol) != dcol: continue
        rv.append(((icol - jcol)[0],(irow - jrow)[0]))
    return rv
def from_npz(filename, pols=None, bls=None, ants=None, verbose=False):
    '''Reconstitute results from to_npz, returns meta, gains, vismdl, xtalk, each
    keyed first by polarization, and then by bl/ant/keyword.
    Optional variables:
    pols: list of polarizations. default: None, return all
    bls: list of baselines. default: None, return all
    ants: list of antennas for gain. default: None, return all
    '''
    if type(filename) is str: filename = [filename]
    if type(pols) is str: pols = [pols]
    if type(bls) is tuple and type(bls[0]) is int: bls = [bls]
    if type(ants) is int: ants = [ants]
    #filename = np.array(filename)
    meta, gains, vismdl, xtalk = {}, {}, {}, {}
    def parse_key(k):
        bl,pol = k.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        return pol,bl
    for f in filename:
        if verbose: print 'Reading', f
        npz = np.load(f)
        for k in npz.files:
            if k[0].isdigit():
                pol,ant = k[-1:],int(k[:-1])
                if (pols==None or pol in pols) and (ants==None or ant in ants): 
                    if not gains.has_key(pol): gains[pol] = {}
                    gains[pol][ant] = gains[pol].get(ant,[]) + [np.copy(npz[k])]
            try: pol,bl = parse_key(k)
            except(ValueError): continue
            if (pols is not None) and (pol not in pols): continue
            if (bls is not None) and (bl not in bls): continue
            if k.startswith('<'):
                if not vismdl.has_key(pol): vismdl[pol] = {}
                vismdl[pol][bl] = vismdl[pol].get(bl,[]) + [np.copy(npz[k])]
            elif k.startswith('('):
                if not xtalk.has_key(pol): xtalk[pol] = {}
                try:
                    dat = np.resize(np.copy(npz[k]),vismdl[pol][vismdl[pol].keys()[0]][0].shape) #resize xtalk to be like vismdl (with a time dimension too)
                except(KeyError):
                    for tempkey in npz.files: 
                        if tempkey.startswith('<'): break
                    dat = np.resize(np.copy(npz[k]),npz[tempkey].shape) #resize xtalk to be like vismdl (with a time dimension too)
                if xtalk[pol].get(bl) is None: #no bl key yet
                    xtalk[pol][bl] = dat
                else: #append to array
                    xtalk[pol][bl] = np.vstack((xtalk[pol].get(bl),dat))
        # for k in [f for f in npz.files if f.startswith('<')]:
        #     pol,bl = parse_key(k)
        #     if not vismdl.has_key(pol): vismdl[pol] = {}
        #     vismdl[pol][bl] = vismdl[pol].get(bl,[]) + [np.copy(npz[k])]
        # for k in [f for f in npz.files if f.startswith('(')]:
        #     pol,bl = parse_key(k)
        #     if not xtalk.has_key(pol): xtalk[pol] = {}
        #     dat = np.resize(np.copy(npz[k]),vismdl[pol][vismdl[pol].keys()[0]][0].shape) #resize xtalk to be like vismdl (with a time dimension too)
        #     if xtalk[pol].get(bl) is None: #no bl key yet
        #         xtalk[pol][bl] = dat
        #     else: #append to array
        #         xtalk[pol][bl] = np.vstack((xtalk[pol].get(bl),dat))
        # for k in [f for f in npz.files if f[0].isdigit()]:
        #     pol,ant = k[-1:],int(k[:-1])
        #     if not gains.has_key(pol): gains[pol] = {}
        #     gains[pol][ant] = gains[pol].get(ant,[]) + [np.copy(npz[k])]
        kws = ['chi','hist','j','l','f']
        for kw in kws:
            for k in [f for f in npz.files if f.startswith(kw)]:
                meta[k] = meta.get(k,[]) + [np.copy(npz[k])]
    #for pol in xtalk: #this is already done above now
        #for bl in xtalk[pol]: xtalk[pol][bl] = np.concatenate(xtalk[pol][bl])
    for pol in vismdl:
        for bl in vismdl[pol]: vismdl[pol][bl] = np.concatenate(vismdl[pol][bl])
    for pol in gains:
        for bl in gains[pol]: gains[pol][bl] = np.concatenate(gains[pol][bl])
    for k in meta:
        try: meta[k] = np.concatenate(meta[k])
        except(ValueError): pass
    return meta, gains, vismdl, xtalk



CONJ = [
    (0,103) , #  1
    (0,95),
    (1,4) ,   #  2
    (1,48),
    (1,18),
    (0,101) , #  3
    (0,62) ,  #  4
    (0,100) , #  5
    (0,97) ,  # 11
    (12,43) , # 12
    (57,64) ] # 15

#SEPS = [(1,4), (1,48),(1,18)]
#SEPS = [(1,4), (1,48)]
SEPS = [(0,103), (0,95)]
#SEPS = [(0,103)]

#CH0,NCHAN = 110, 51
CH0,NCHAN = 125, 21
#bandpass = np.load('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/bandpass.npz')['bandpass']
bandpass = np.load('bandpass.npz')['bandpass']
bandpass.shape = (1,-1)

cwd = os.getcwd()
if cwd.startswith('/Users/yunfanzhang/'):
    dataDIR = '/Users/yunfanzhang/local/DATA128/DATA/'
elif cwd.startswith('/Users/yunfanz/'):
    dataDIR = '/Users/yunfanz/Data/PAPER128/DATA/'
elif cwd.startswith('/home/yunfanz/'):
    #dataDIR = '/data2/PAPER/omni_v2_xtalk/'
    dataDIR = '/data2/PAPER/2013_e1_omv3/'

sets = {}
for fn in np.sort(glob.glob(dataDIR+'zen*.xx.npz')):
    k = fn.split('.')[1]
    sets[k] = sets.get(k,[]) + [fn]

data,wgts = {}, {}
lsts = {}
chisqs = {}
for s in sets:
    if not lsts.has_key(s):
        meta, gains, vismdl, xtalk = from_npz(sets[s], bls=SEPS, pols='xx', ants=1,verbose=True)
        lsts[s] = meta['lsts']
    chisqs[s] = meta['chisq'][:,CH0:CH0+NCHAN]
    for pol in vismdl:
        for bl in SEPS:
            k = (s,pol,bl)
            data[k] = vismdl[pol][bl][:,CH0:CH0+NCHAN]
            if bl in CONJ: data[k] = data[k].conj()
            data[k] *= bandpass[:,CH0:CH0+NCHAN]
            wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1)

ind = {}

#import IPython; IPython.embed()
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#inds = oqe.lst_align(lsts, lstres=lst_res)
#data,wgts = oqe.lst_align_data(inds, dsets=data, wgts=wgts)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# for s in sets: chisqs[s] = chisqs[s][ind[s]].T
########################################################################
    

ks = [(s,'xx',bl) for bl in SEPS for s in sets]
NK = len(ks)


data_gt, wgt_gt = {},{}
nlst = 2400
for k in data:
    lst_g,data_gt[k],wgt_gt[k] = oqe.lst_grid(lsts[k[0]],data[k],lstbins=nlst)
    data_gt[k], wgt_gt[k] = data_gt[k][nlst/3:nlst*4/5], wgt_gt[k][nlst/3:nlst*4/5]
    #wgt_gt[k] = np.where(wgt_gt[k]>0.5*np.max(wgt_gt[k]),1,0)
    lst_g = lst_g[nlst/3:nlst*4/5]

set1,set2 = 'even', 'odd'
def is_odd(num):
    return num & 0x1

data_g, wgt_g = {},{}
for bl in SEPS:
    k1 = ('odd','xx',bl)
    k2 = ('even','xx',bl)

    for k in data_gt.keys():
        if is_odd(int(k[0])) and k[1]==k1[1] and k[2]==k1[2]:  #i.e. k satisfies k1
            temp = data_gt[k]; temp.shape += (1,)
            data_g[k1] = data_g.get(k1,[])+[temp]
            wtemp = wgt_gt[k]; wtemp.shape += (1,)
            wgt_g[k1] = wgt_g.get(k1,[])+[wtemp]
        elif (not is_odd(int(k[0]))) and k[1]==k2[1] and k[2]==k2[2]:
            temp = data_gt[k]; temp.shape += (1,)
            data_g[k2] = data_g.get(k2,[])+[temp]
            wtemp = wgt_gt[k]; wtemp.shape += (1,)
            wgt_g[k2] = wgt_g.get(k2,[])+[wtemp]
    #import IPython; IPython.embed()
    for k in (k1, k2):
        data_g[k] = np.dstack(data_g[k])
        wgt_g[k] = np.dstack(wgt_g[k])
        data_g[k] = np.ma.average(data_g[k], axis=2, weights=wgt_g[k])
        data_g[k] = data_g[k].filled(fill_value=0)
        wgt_g[k] = np.mean(wgt_g[k], axis=2)
        wgt_g[k] = np.where(wgt_g[k]>0.5*np.max(wgt_g[k]),1,0)
    #import IPython; IPython.embed()

    print k1, k2
    print 'should be even_xx_0_95 odd_xx_0_103'
# np.savez('griddata14148', even_xx_1_4=data_g[k1], odd_xx_1_48=data_g[k2])
# np.savez('gridwgt14148', even_xx_1_4=wgt_g[k1], odd_xx_1_48=wgt_g[k2])
# np.savez('gridlst14148', lst=lst_g)
np.savez('griddata2013', odd_xx_0_103=data_g[('odd','xx',(0,103))], even_xx_0_103=data_g[('even','xx',(0,103))],odd_xx_0_95=data_g[('odd','xx',(0,95))], even_xx_0_95=data_g[('even','xx',(0,95))])
np.savez('gridwgt2013', odd_xx_0_103=wgt_g[('odd','xx',(0,103))], even_xx_0_103=wgt_g[('even','xx',(0,103))],odd_xx_0_95=wgt_g[('odd','xx',(0,95))], even_xx_0_95=wgt_g[('even','xx',(0,95))])
np.savez('gridlst2013', lst=lst_g)
import IPython; IPython.embed()
