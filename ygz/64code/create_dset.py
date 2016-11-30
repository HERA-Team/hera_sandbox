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

SEPS = [(1,4), (1,48),(1,18)]
SEPS = [(1,4), (1,48)]
SEPS = [(0,103), (0,95)]

# CH0,NCHAN = 110, 51
# #bandpass = np.load('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/bandpass.npz')['bandpass']
# bandpass = np.load('bandpass.npz')['bandpass']
# bandpass.shape = (1,-1)
# fqs = np.linspace(.1,.2,bandpass.size)
# WINDOW = 'blackman-harris'
# win = aipy.dsp.gen_window(NCHAN, WINDOW); win.shape = (-1,1)
# afreqs = fqs[CH0:CH0+NCHAN]
# fq = np.average(afreqs)
# z = capo.pspec.f2z(fq)
# sdf = fqs[1]-fqs[0]
# B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW['none'] #proper normalization
# etas = np.fft.fftshift(capo.pspec.f2eta(afreqs)) #create etas (fourier dual to frequency)
# kpl = etas * capo.pspec.dk_deta(z)
# bm = np.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq) * 2.35 # correction for beam^2
# scalar = capo.pspec.X2Y(z) * bm * B
# B_win = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] #proper normalization
# scalar_win = capo.pspec.X2Y(z) * bm * B_win

#dataDIR = '/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/'
cwd = os.getcwd()
if cwd.startswith('/Users/yunfanzhang/'):
    dataDIR = '/Users/yunfanzhang/local/DATA128/DATA/'
elif cwd.startswith('/Users/yunfanz/'):
    dataDIR = '/Users/yunfanz/Data/PAPER128/DATA/'
elif cwd.startswith('/home/yunfanz/'):
    dataDIR = '/home/yunfanz/Projects/21cm/Data/DATA128/DATA/'

sets = {}
for fn in glob.glob(dataDIR+'zen*.xx.npz'):
    k = fn.split('.')[1]
    sets[k] = sets.get(k,[]) + [fn]

# sets = {
#     #'day0' : sys.argv[1:],
#     'day0' : glob.glob(dataDIR+'zen.2456710.*.xx.npz'),
#     'day1' : glob.glob(dataDIR+'zen.2456715.*.xx.npz'),
#     'day2' : glob.glob(dataDIR+'zen.2456716.*.xx.npz'),
# }
data,wgts = {}, {}
lsts = {}
chisqs = {}
for s in sets:
    if not lsts.has_key(s):
        meta, gains, vismdl, xtalk = capo.omni.from_npz(sets[s], bls=SEPS, pols='xx', ants=1,verbose=True)
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
    # lst_g,data_g[k],wgt_g[k] = oqe.lst_grid(lsts[k[0]],data[k],lstbins=6000)
    # data_g[k], wgt_g[k] = data_g[k][2200:5000], wgt_g[k][2200:5000]
    lst_g,data_gt[k],wgt_gt[k] = oqe.lst_grid(lsts[k[0]],data[k],lstbins=nlst)
    data_gt[k], wgt_gt[k] = data_gt[k][nlst/3:nlst*4/5], wgt_gt[k][nlst/3:nlst*4/5]
    wgt_gt[k] = np.where(wgt_gt[k]>0.5*np.max(wgt_gt[k]),1,0)
# dlst = np.average(lst_g[1:] - lst_g[:-1])/2
#import IPython; IPython.embed()
set1,set2 = 'even', 'odd'
def is_odd(num):
    return num & 0x1
k1 = ('odd','xx',(0,103))
k2 = ('even','xx',(0,95))
data_g, wgt_g = {},{}
for k in data_gt.keys():
    if is_odd(int(k[0])) and k[1]==k1[1] and k[2]==k1[2]:
        temp = data_gt[k]; temp.shape += (1,)
        data_g[k1] = data_g.get(k1,[])+[temp]
        wtemp = wgt_gt[k]; wtemp.shape += (1,)
        wgt_g[k1] = wgt_g.get(k1,[])+[wtemp]
    elif (not is_odd(int(k[0]))) and k[1]==k2[1] and k[2]==k2[2]:
        temp = data_gt[k]; temp.shape += (1,)
        data_g[k2] = data_g.get(k2,[])+[temp]
        wtemp = wgt_gt[k]; wtemp.shape += (1,)
        wgt_g[k2] = wgt_g.get(k2,[])+[wtemp]
data_g[k1] = np.mean(np.dstack(data_g[k1]), axis=2)
data_g[k2] = np.mean(np.dstack(data_g[k2]), axis=2)
wgt_g[k1] = np.mean(np.dstack(wgt_g[k1]), axis=2)
wgt_g[k2] = np.mean(np.dstack(wgt_g[k2]), axis=2)
print k1, k2
print 'should be even_xx_0_95 odd_xx_0_103'
#np.savez('griddata', even_xx_0_95=data_g[k1], odd_xx_0_103=data_g[k2])
#np.savez('gridwgt', even_xx_0_95=wgt_g[k1], odd_xx_0_103=wgt_g[k2])
np.savez('gridlst', lst=lst_g)
import IPython; IPython.embed()
