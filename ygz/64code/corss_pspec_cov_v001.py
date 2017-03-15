#! /usr/bin/env python

import numpy as np, aipy, capo, pylab as plt, sys, glob
import md5
from joblib import Parallel, delayed
import multiprocessing, joblib
import oqe
import os
from waterfall import waterfall
num_cores = multiprocessing.cpu_count()

"""See create_dset"""
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

# ks = [(s,'xx',bl) for bl in SEPS for s in sets]
# NK = len(ks)
def parse_string(s):
    p = s.split('_')
    k = (p[0], p[1], (int(p[2]), int(p[3])))
    return k

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

def get_p(k1,k2,mode,offset=0, rephs=0):

    assert(mode in 'IWC')
    if True:
            #save_data,save_wgt = {},{}
        ds_new = None
        data_g[k1] *= np.exp(1j*rephs*afreqs)[np.newaxis, :]
        #import IPython; IPython.embed()
        if offset >1e-8:
            ds_new = oqe.DataSet({k1:data_g[k1][:-offset], k2:data_g[k2][offset:]}, {k1:wgt_g[k1][:-offset], k2:wgt_g[k2][offset:]})
        elif offset < -1e-8:
            ds_new = oqe.DataSet({k1:data_g[k1][-offset:], k2:data_g[k2][:offset]}, {k1:wgt_g[k1][-offset:], k2:wgt_g[k2][:offset]})
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
        if offset > 0:
            # k1:   (#####)##
            # k2: ##(#####)
            of_data1 = data_g[k1][:-offset].T; of_data2 = data_g[k2][offset:].T
        elif offset < 0:
            of_data1 = data_g[k1][-offset:].T; of_data2 = data_g[k2][:offset].T
        else:
            of_data1 = data_g[k1].T; of_data2 = data_g[k2].T
        
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
offset_dict[((0,103), (0,103))] = 0.0
rephs_dict = {((0,103), (0,103)): 0.0, ((0,103),(0,95)): 104.936, ((0,95),(0,103)): -104.936}
#ind[set1], ind[set2] = lst_align(lsts[set1], lsts[set2])

num = '2013'
from itertools import product
d_file = np.load('griddata'+num+'.npz')
w_file = np.load('gridwgt'+num+'.npz')
data_g, wgt_g, K = {},{}, []
for s in d_file.files:
    k = parse_string(s)
    K.append(k)
    data_g[k] = d_file[s]
    wgt_g[k] = w_file[s]
lst_g = np.load('gridlst'+num+'.npz')['lst']
dlst = lst_g[1]-lst_g[0]
################################
ds = oqe.DataSet(data_g, wgt_g)
#import IPython; IPython.embed()
set_C(ds,3e4)
#set_C(0)
#pI,pW,pC = get_p(ks[0],ks[1])
#################################

oflist = (np.arange(400)-200)#+int(0.0548/dlst)
oflst = oflist*dlst
#oflist = oflist[oflist>=0]
n_jobs = 8
# of_batches = []
# for n in xrange(n_jobs):
#     of_batches.append(oflist[n::n_jobs])
print K
k1,k2 = K[0], K[3]
print k1, k2
rephs = rephs_dict[(k1[-1], k2[-1])]
#import IPython; IPython.embed()

# d = data_g[k1].T
# d2 = data_g[k2].T
# d_ = np.fft.fft(d, axis=0); d2_ = np.fft.fft(d2, axis=0)
# Cd = np.fft.ifftshift(d2_.conj()*d_,axes=0)
# plt.figure()
# plt.imshow(np.abs(Cd.T), aspect='auto')
# plt.show()
# import IPython; IPython.embed()




res = Parallel(n_jobs=n_jobs)(delayed(get_p)(k1,k2,'I',ofbatch, rephs=rephs) for ofbatch in oflist)
# def postprocess(res):
#     pC, OFST = zip(*res)
#     OFST = np.hstack(OFST)
#     pC = np.hstack(pC)
#     inds = sorted(range(len(OFST)), key=lambda k: OFST[k])
#     return OFST[inds], pC[inds]
# #ofst, pC = postprocess(res)
pC, OFST = zip(*res)
inds = sorted(range(len(OFST)), key=lambda k: pC[k].shape[1])
Nlst = pC[inds[0]].shape[1]
#import IPython; IPython.embed()
pCt = np.zeros((oflist.size, NCHAN, Nlst), dtype=pC[0].dtype)

for i, psc in enumerate(pC):
    pCt[i] = psc[:,-Nlst:]
ps = np.mean(pCt, axis=2)
#import IPython; IPython.embed()
plt.figure()
waterfall(ps.T, mx=16, drng=7)
plt.colorbar()
plt.show()
# for cnt, ofst in enumerate(OFST):
#     psc = pC[cnt]
#     plt.subplot(10,1,cnt+1)
#     plt.title(ofst)
#     capo.plot.waterfall(psc, mx=16, drng=7)
#     plt.colorbar()
# plt.show()
import IPython; IPython.embed()
