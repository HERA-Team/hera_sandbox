#! /usr/bin/env python

import numpy as np, aipy, capo, pylab as plt, sys, glob
import md5
from joblib import Parallel, delayed
import multiprocessing
import oqe
import os
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

def rebin_lst(binsize, lsts, d, w):  #d: data in time by freq. 
    bins = lsts/binsize
    b0 = int(bins[0])
    bmax = int(np.ceil(np.max(bins)))
    dbin = np.zeros((bmax,d.shape[1]), dtype=d.dtype)
    dwgt = np.zeros(dbin.shape, dtype=np.float)
    for i,b in enumerate(bins-b0):
        db = b % 1
        w1,w2 = (1-db) * w[i], db * w[i] # linearly interpolate between bins
        dbin[np.floor(b)] += w1 * d[i]
        dwgt[np.floor(b)] += w1
        dbin[np.ceil(b)] += w2 * d[i]
        dwgt[np.ceil(b)] += w2
    lstbin = (b0 + np.arange(bmax)) * binsize
    return lstbin, dbin / np.where(dwgt > 0, dwgt, 1), dwgt


CONJ = [
    (0,103) , #  1
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
#SEPS += [(2,105), (1,83)]
#SEPS += [(0,79), (0,78)]
#SEPS += [(0,70),(0,71)]
#SEPS += [(1,107),(0,51)]
#SEPS += [(3,105),(3,106)]
#CH0,NCHAN = 90, 31
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
elif cwd.startswith('/home/yunfanz/'):
    dataDIR = '/home/yunfanz/EoR/DATA128/DATA/'
sets = {
    #'day0' : sys.argv[1:],
    #'day0' : glob.glob('zen.2456714.*.xx.npz'),
    'day1' : glob.glob(dataDIR+'zen.2456715.5*.xx.npz'),
    'day2' : glob.glob(dataDIR+'zen.2456716.5*.xx.npz'),
}
data,wgts = {}, {}
lsts = {}
chisqs = {}
# def from_npz(file):
#      res = capo.omni.from_npz(file,verbose=True)
#      return res
for s in sets:
    if not lsts.has_key(s):
        # res = Parallel(n_jobs=4)(delayed(from_npz)(file) for file in sets[s])
        # #import IPython; IPython.embed()
        # meta, gains, vismdl, xtalk = [],[],[],[]
        # for i, elt in enumerate(res):
        #     meta.append(elt[0])
        #     gains.append(elt[1])
        #     vismdl.append(elt[2])
        #     xtalk.append(elt[3])
        meta, gains, vismdl, xtalk = capo.omni.from_npz(sets[s], bls=SEPS, pols='xx', ants=1,verbose=True)
        lsts[s] = meta['lsts']
    chisqs[s] = meta['chisq'][:,CH0:CH0+NCHAN]
    for pol in vismdl:
        #for bl in vismdl[pol]:
        for bl in SEPS:
            k = (s,pol,bl)
            data[k] = vismdl[pol][bl][:,CH0:CH0+NCHAN]
            if bl in CONJ: data[k] = data[k].conj()
            data[k] *= bandpass[:,CH0:CH0+NCHAN]
            wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1)

ind = {}
set1,set2 = sets.keys()[0], sets.keys()[-1]
lst_res = np.average(lsts[set1][1:] - lsts[set1][:-1])/2
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#inds = oqe.lst_align(lsts, lstres=lst_res)
#data,wgts = oqe.lst_align_data(inds, dsets=data, wgts=wgts)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# for s in sets: chisqs[s] = chisqs[s][ind[s]].T
########################################################################
    

ks = [(s,'xx',bl) for bl in SEPS for s in sets]
#k1a, = [(s,'xx',(0,103)) for s in sets]
#k2a, = [(s,'xx',(0,111)) for s in sets]
#k3a, = [(s,'xx',(0, 95)) for s in sets]
#ks = [k1a,k2a,k3a]
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
        #import IPython; IPython.embed()
        #w = np.where(ds.w[ki] > 0, 1, 0)
        #Cs[k] = sum([oqe.cov(ds.x[ki][:,400:],w[:,400:])+norm*np.identity(NCHAN) for ki in ks if ki[2] != k[2]])
        #Cs[k] = sum([ds.C(k)+norm*np.identity(NCHAN) for ki in ks if ki != k])
        #Cs[k] = sum([ds.C(k)+norm*np.identity(NCHAN) for ki in data if ki[2] != k[2]])
        dst.set_C({k:Cs[k]})
        iCs[k] = dst.iC(k)

def get_p(k1,k2,mode,offset=0):

    assert(mode in 'IWC')
    if mode == 'I':
        qI = ds.q_hat(k1,k2,use_cov=False)
        FI = ds.get_F(k1,k2,use_cov=False)
        MI,WI = ds.get_MW(FI, mode='I')
        pI = ds.p_hat(MI,qI)
        return pI * scalar
    elif mode == 'W':
        pW = 1.6*2*np.fft.fftshift(np.fft.ifft(win*data[k1].T, axis=0) * np.fft.ifft(win*data[k2].T, axis=0).conj(), axes=0)
        return pW * scalar_win
    elif mode == 'C':
        if True:
            #save_data,save_wgt = {},{}
            ds_new = None
            if abs(offset) >1e-8:
                # save_data,save_wgt = data_g.copy(),wgt_g.copy()
                # #save_data,save_wgt ={},{}
                # save_data[k1] = data_g[k1][ind_offset:]; save_wgt[k1] = wgt_g[k1][ind_offset:]
                # if ind_offset>0: save_data[k2] = data_g[k2][:-ind_offset]; save_wgt[k2] = wgt_g[k2][:-ind_offset]
                # ds_new = oqe.DataSet(save_data, save_wgt)
                #import IPython; IPython.embed()

                ds_new = oqe.DataSet({k1:data_g[k1][:-offset], k2:data_g[k2][offset:]}, {k1:wgt_g[k1][:-offset], k2:wgt_g[k2][offset:]})
            else:
                ds_new = ds
            # iC_dict = {}
            # for k in ds.x:
            #     iC_dict[k] = ds.iC(k)
            # ds_new.set_iC(iC_dict)
            ds_new.set_C(ds._C)
            #import IPython; IPython.embed()
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
            print "computing power spectrum for ", k1, k2
            qC = ds_new.q_hat(k1,k2)
            FC = ds_new.get_F(k1,k2)
            #import IPython; IPython.embed()
            MC,WC = ds_new.get_MW(FC, mode='F^-1/2')
            pC = ds_new.p_hat(MC,qC)
        return pC * scalar

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
dlst = lst_res
#ind[set1], ind[set2] = lst_align(lsts[set1], lsts[set2])

from itertools import product
# for bl1, bl2 in product(SEPS,SEPS):
#     if bl1 == bl2: offset = 0
#     else: offset = offset_dict[(bl1,bl2)]
#     ind_offset = int(offset/dlst)
#     for pol in vismdl:
#         for s1,s2 in product(sets,sets):
#             if s1 == s2: continue


# for k in data.keys():
#     (s,pol,bl) = k
#     data[k] = data[k][ind[s]]     #aligning the lsts
#     wgts[k] = wgts[k][ind[s]]

data_g, wgt_g = {},{}
for k in data:
    lst_g,data_g[k],wgt_g[k] = oqe.lst_grid(lsts[k[0]],data[k])
################################
ds = oqe.DataSet(data_g, wgt_g)
#import IPython; IPython.embed()
set_C(ds,3e-6)
#set_C(0)
#pI,pW,pC = get_p(ks[0],ks[1])
#################################
for cnt,k in enumerate(ks):
    plt.subplot(NK,1,cnt+1)
    capo.plot.waterfall(ds.x[k], drng=3)
    plt.colorbar()
#plt.savefig('fig1.png')
#plt.show()


sep_pairs = product(SEPS,SEPS)
for cnt, bls in enumerate(sep_pairs):
    k1 = (set1,pol,bls[0])
    k2 = (set2,pol,bls[1])
    if bls[0] == bls[1]: offset = 0
    else: offset = offset_dict[(bls[0],bls[1])]
    ind_offset = int(offset/dlst)
    print ind_offset
    pC = get_p(k1,k2,'C',offset=ind_offset)
    plt.subplot(5,1,cnt+1)
    plt.title(bls)
    capo.plot.waterfall(pC, mx=16, drng=7)
    plt.colorbar()
plt.show()
import IPython; IPython.embed()
#for cnt,k in enumerate(ks):
#    print k
    #plt.subplot(NK,1,cnt+1)
    #pC = get_p(k,k,'C')
    #plt.title(k[0])
    #capo.plot.waterfall(pC, mx=16, drng=7)
    #plt.colorbar()
#plt.savefig('fig2.png')
#plt.show()

'''
pC1 = get_p(k1a,k1b,'C')
pC2 = get_p(k2a,k2b,'C')
pC3 = get_p(k3a,k3b,'C')
plt.plot(np.abs(np.median(pC1, axis=1).real), 'r', label='pC1',)
plt.plot(np.abs(np.median(pC2, axis=1).real), 'b', label='pC2',)
plt.plot(np.abs(np.median(pC3, axis=1).real), 'k', label='pC3',)
plt.legend()
plt.show()
'''
