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

CH0,NCHAN = 110, 51
#bandpass = np.load('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/bandpass.npz')['bandpass']
bandpass = np.load('bandpass.npz')['bandpass']
bandpass.shape = (1,-1)
fqs = np.linspace(.1,.2,bandpass.size)
WINDOW = 'blackman-harris'
win = aipy.dsp.gen_window(NCHAN, WINDOW); win.shape = (-1,1)
afreqs = fqs[CH0:CH0+NCHAN]
fq = np.average(afreqs)
z = capo.pspec.f2z(fq)
sdf = fqs[1]-fqs[0]
B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW['none'] #proper normalization
etas = np.fft.fftshift(capo.pspec.f2eta(afreqs)) #create etas (fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z)
bm = np.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq) * 2.35 # correction for beam^2
scalar = capo.pspec.X2Y(z) * bm * B
B_win = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] #proper normalization
scalar_win = capo.pspec.X2Y(z) * bm * B_win

#dataDIR = '/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/'
cwd = os.getcwd()
if cwd.startswith('/Users/yunfanzhang/'):
    dataDIR = '/Users/yunfanzhang/local/DATA128/DATA/'
elif cwd.startswith('/Users/yunfanz/'):
    dataDIR = '/Users/yunfanz/Data/PAPER128/DATA/'
elif cwd.startswith('/home/yunfanz/'):
    dataDIR = '/home/yunfanz/Projects/21cm/Data/DATA128/DATA/'
sets = {
    #'day0' : sys.argv[1:],
    #'day0' : glob.glob('zen.2456714.*.xx.npz'),
    'day1' : glob.glob(dataDIR+'zen.2456715.*.xx.npz'),
    'day2' : glob.glob(dataDIR+'zen.2456716.*.xx.npz'),
}
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
set1,set2 = sets.keys()[0], sets.keys()[-1]

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#inds = oqe.lst_align(lsts, lstres=lst_res)
#data,wgts = oqe.lst_align_data(inds, dsets=data, wgts=wgts)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# for s in sets: chisqs[s] = chisqs[s][ind[s]].T
########################################################################
    

ks = [(s,'xx',bl) for bl in SEPS for s in sets]
NK = len(ks)

def set_C(dst,norm=3e-6):
    dst.clear_cache()
    Cs,iCs = {},{}
    #import IPython; IPython.embed()
    for k in dst.x.keys():
        #Cs[k] = sum([oqe.cov(ds.x[k][:,400:],ds.w[k][:,400:])+norm*np.identity(NCHAN) for ki in ks if ki != k])
        #Cs[k] = sum([oqe.cov(ds.x[ki][:,400:],ds.w[ki][:,400:])+norm*np.identity(NCHAN) for ki in ks if ki[2] != k[2]])
        Ndim = dst.x[k].shape[0]
        Cs[k] = oqe.cov(dst.x[k][:,0:],dst.w[k][:,0:])+norm*np.identity(Ndim)
        dst.set_C({k:Cs[k]})
        iCs[k] = dst.iC(k)

def get_p(k1,k2,mode,offset=0):

    assert(mode in 'IWC')
    if True:
            #save_data,save_wgt = {},{}
        ds_new = None
        if abs(offset) >1e-8:
            ds_new = oqe.DataSet({k1:data_g[k1][:-offset], k2:data_g[k2][offset:]}, {k1:wgt_g[k1][:-offset], k2:wgt_g[k2][offset:]})
        else:
            ds_new = ds
        ds_new.set_C(ds._C)
        #import IPython; IPython.embed()
    if mode == 'I':
        qI = ds_new.q_hat(k1,k2,use_cov=False)
        FI = ds_new.get_F(k1,k2,use_cov=False)
        MI,WI = ds_new.get_MW(FI, mode='I')
        pI = ds_new.p_hat(MI,qI)
        return pI * scalar, offset
    elif mode == 'W':
        of_data1 = data_g[k1][:-offset].T; of_data2 = data_g[k2][offset:].T
        pW = 1.6*2*np.fft.fftshift(np.fft.ifft(win*of_data1, axis=0) * np.fft.ifft(win*of_data2, axis=0).conj(), axes=0)
        return pW * scalar_win, offset
    elif mode == 'C':
        if False:
            save_iC = {}
            for k in (k1,k2): save_iC[k] = ds.iC(k).copy()
            iCs,sums = {},{}
            for k in (k1,k2):
                iCs[k] = {}
                sums[k] = [md5.md5(ds.w[k][:,i]).digest() for i in xrange(ds.w[k].shape[1])]
                for i,s in enumerate(sums[k]): iCs[k][s] = ds.w[k][:,i:i+1]
                for s in iCs[k].keys():
                    w = iCs[k][s]
                    iCs[k][s] = np.linalg.pinv(ds.C(k) * (w * w.T), rcond=1e-12)
            pCs = {}
            for s in iCs[k].keys():
                ds.set_iC({k1:iCs[k1][s], k2:iCs[k2][s]})
                qC = ds.q_hat(k1,k2)
                FC = ds.get_F(k1,k2)
                MC,WC = ds.get_MW(FC, mode='F^-1/2')
                pCs[s] = ds.p_hat(MC,qC)
            ds.set_iC(save_iC)
            # XXX deal with diff w for k1,k2
            pC = np.array([pCs[sums[k1][i]][:,i] for i in xrange(ds.w[k1].shape[1])]).T
        else:
            print "computing power spectrum for ", k1, k2, offset
            qC = ds_new.q_hat(k1,k2)
            FC = ds_new.get_F(k1,k2)
            #import IPython; IPython.embed()
            MC,WC = ds_new.get_MW(FC, mode='F^-1/2')
            pC = ds_new.p_hat(MC,qC)
        return pC * scalar, offset
def get_p_batch(k1,k2,mode,offset_ls=[0]):
    respC = []
    for of in offset_ls:
        pc, _ = get_p(k1,k2,mode,of)
    respC.append(pc)
    return respC, offset_ls

def lst_align(lsts1,lsts2,offset=0):
    i=0
    if lsts1[0]+offset>lsts2[0]:
        while lsts1[0]+offset>lsts2[i]: i += 1
        Npt = min(lsts1.size,lsts2.size-i)
        return np.arange(0,Npt),np.arange(i,Npt+i)
    else:
        while lsts1[i]+offset>lsts2[0]: i += 1
        Npt = min(lsts1.size-i,lsts2.size)
        return np.arange(i,Npt+i),n.arange(0,Npt)

offset_dict = {((1,48), (1,4)):0.031,((1,4), (1,48)):0.031}
offset_dict[((0,103),(0,95))] = 0.0548
offset_dict[((0,95),(0,103))] = 0.0548

#ind[set1], ind[set2] = lst_align(lsts[set1], lsts[set2])

from itertools import product

data_g, wgt_g = {},{}
nlst = 2400
for k in data:
    # lst_g,data_g[k],wgt_g[k] = oqe.lst_grid(lsts[k[0]],data[k],lstbins=6000)
    # data_g[k], wgt_g[k] = data_g[k][2200:5000], wgt_g[k][2200:5000]
    lst_g,data_g[k],wgt_g[k] = oqe.lst_grid(lsts[k[0]],data[k],lstbins=nlst)
    data_g[k], wgt_g[k] = data_g[k][nlst/3:nlst*4/5], wgt_g[k][nlst/3:nlst*4/5]
    wgt_g[k] = np.where(wgt_g[k]>0.5*np.max(wgt_g[k]),1,0)
dlst = np.average(lst_g[1:] - lst_g[:-1])/2
################################
ds = oqe.DataSet(data_g, wgt_g)
#import IPython; IPython.embed()
set_C(ds,3e4)
#set_C(0)
#pI,pW,pC = get_p(ks[0],ks[1])
#################################


k1 = (set1,'xx',(0,103))
k2 = (set1,'xx',(0,95))
oflist = (np.arange(100)-50)+int(0.0548/dlst)
oflist = oflist[oflist>=0]
n_jobs = num_cores
# of_batches = []
# for n in xrange(n_jobs):
#     of_batches.append(oflist[n::n_jobs])
res = Parallel(n_jobs=n_jobs)(delayed(get_p)(k1,k2,'C',ofbatch) for ofbatch in oflist)
# def postprocess(res):
#     pC, OFST = zip(*res)
#     OFST = np.hstack(OFST)
#     pC = np.hstack(pC)
#     inds = sorted(range(len(OFST)), key=lambda k: OFST[k])
#     return OFST[inds], pC[inds]
# #ofst, pC = postprocess(res)
pC, OFST = zip(*res)
Nlst = pC[-1].shape[1]
pCt = np.zeros((oflist.size, NCHAN, Nlst))
import IPython; IPython.embed()
for i, psc in enumerate(pC):
    pCt[i] = psc[:,-Nlst:]
#waterfall(np.mean(pCt, axis=2), mx=16, drng=7)

# for cnt, ofst in enumerate(OFST):
#     psc = pC[cnt]
#     plt.subplot(10,1,cnt+1)
#     plt.title(ofst)
#     capo.plot.waterfall(psc, mx=16, drng=7)
#     plt.colorbar()
# plt.show()
#import IPython; IPython.embed()
