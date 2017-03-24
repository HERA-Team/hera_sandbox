#! /usr/bin/env python
# import matplotlib
# matplotlib.use('TkAgg')
import numpy as np, aipy, capo, matplotlib.pyplot as plt, sys, glob
import md5
import oqe, os
import fileinput
import cProfile
from waterfall import waterfall

pr = cProfile.Profile()
def dB(sig): return 10*np.log10(np.abs(np.average(sig.real, axis=1)))

def find_sep(aa, bls, drow=None, dcol=None):
    layout = aa.ant_layout
    rv = []
    for i,j in bls:
        irow,icol = np.where(layout == i)
        jrow,jcol = np.where(layout == j)
        if not drow is None and abs(irow - jrow) != drow: continue
        if not dcol is None and abs(icol - jcol) != dcol: continue
        rv.append((i,j))
    return rv

def rebin_lst(binsize, lsts, d, w):
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

SEPS = [
    (0,103) , #  1
    (1,4) ,   #  2
    (0,101) , #  3
    (0,62) ,  #  4
    (0,100) , #  5
    (1,13) ,  #  6
    (1,70) ,  #  7
    (1,56) ,  #  8
    (1,71) ,  #  9
    (1,59) ,  # 10
    (0,97) ,  # 11
    (12,43) , # 12
    (9,71) ,  # 13
    (9,59) ,  # 14
    (57,64) , # 15
]

CONJ = [
    (0,103) , #  1
    (1,4) ,   #  2
    (0,101) , #  3
    (0,62) ,  #  4
    (0,100) , #  5
    (0,97) ,  # 11
    (12,43) , # 12
    (57,64) ] # 15
SEPS = [(0,103)]
CH0,NCHAN = 110, 51
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
cwd = os.getcwd()
if cwd.startswith('/Users/yunfanzhang/'):
    dataDIR = '/Users/yunfanzhang/local/DATA128/DATA/'
elif cwd.startswith('/Users/yunfanz/'):
    dataDIR = '/Users/yunfanz/Data/PAPER128/DATA/'
elif cwd.startswith('/home/yunfanz/'):
    dataDIR = '/data2/PAPER/omni_v2_xtalk/'
sets = {
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
lst_res = np.average(lsts[set1][1:] - lsts[set1][:-1])/2

ks = [(s,'xx',bl) for bl in SEPS for s in sets]

NK = len(ks)

def set_C(norm=30.):
    ds.clear_cache()
    Cs,iCs = {},{}
    for k in ks:
        Ndim = ds.x[k].shape[0]
        Cs[k] = oqe.cov(ds.x[k][:,:],ds.w[k][:,])+norm*np.identity(Ndim)
        ds.set_C({k:Cs[k]})
        iCs[k] = ds.iC(k)

def get_p(k1,k2,mode):
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
            #pr.enable()
            qC = ds.q_hat(k1,k2,cov_flagging=True)
            FC = ds.get_F(k1,k2,cov_flagging=True)
            MC,WC = ds.get_MW(FC, mode='F^-1/2')
            pC = ds.p_hat(MC,qC)
            #pr.disable()
            #pr.print_stats(sort='time')
        return pC * scalar

data_g, wgt_g = {},{}
for k in data:
    lst_g,data_g[k],wgt_g[k] = oqe.lst_grid(lsts[k[0]],data[k],lstbins=1500)
    data_g[k], wgt_g[k] = data_g[k][550:1200], wgt_g[k][550:1200]
    # lst_g,data_g[k],wgt_g[k] = oqe.lst_grid(lsts[k[0]],data[k],lstbins=6000)
    # data_g[k], wgt_g[k] = data_g[k][2200:5000], wgt_g[k][2200:5000]
    wgt_g[k] = np.where(wgt_g[k]>0.5*np.max(wgt_g[k]),1,0)
k1, k2 = data.keys()
data_g[('mean','xx',(0,103))] = (data_g[k1]+data_g[k2])/2
wgt_g[('mean','xx',(0,103))] = (wgt_g[k1]+wgt_g[k2])/2
################################
#import IPython; IPython.embed()

#import IPython; IPython.embed()
# def lst_align(lsts1,lsts2,lstres,offset=0):
#     i=0
#     if lsts1[0]+offset>lsts2[0]:
#         while lsts1[0]+offset>lsts2[i]: i += 1
#         Npt = min(lsts1.size,lsts2.size-i)
#         return np.arange(0,Npt),np.arange(i,Npt+i)
#     else:
#         while lsts1[i]+offset>lsts2[0]: i += 1
#         Npt = min(lsts1.size-i,lsts2.size)
#         return np.arange(i,Npt+i),n.arange(0,Npt)

# ind[set1], ind[set2] = lst_align(lsts[set1], lsts[set2], lstres=lst_res)

# for k in data.keys():
#     (s,pol,bl) = k
#     data[k] = data[k][ind[s]]     #aligning the lsts
#     wgts[k] = wgts[k][ind[s]]
# ds = capo.oqe.DataSet(data, wgts)



ds = oqe.DataSet(data_g, wgt_g)
set_C(norm=30000.)
#import IPython; IPython.embed()
# for cnt,k in enumerate(ks):
#     plt.subplot(NK,1,cnt+1)
#     capo.plot.waterfall(ds.x[k], drng=3)
#     plt.title(k)
#     plt.colorbar()
# plt.savefig('timeseries.png')


bls = (0,103)
k1 = (set1,pol,(0,103))
k2 = (set2,pol,(0,103))
k3 = ('mean','xx',(0,103))

f, (ax1, ax2, ax3) = plt.subplots(3,1)
pC1 = get_p(k1,k1,'C')
plt.title(set2+set2+str(bls)+'C')
im1 = waterfall(pC1, ax=ax1, mx=16, drng=7)
#plt.colorbar()
#plt.figure()
pC2 = get_p(k1,k2,'C')
plt.title(set1+set1+str(bls)+'I')
im2 = waterfall(pC2, ax=ax2, mx=16, drng=7)
#plt.colorbar()
#pC = get_p(k1,k2,'C')
#plt.subplot(3,1,3)
#plt.title(set1+set2+str(bls)+'C')
#capo.plot.waterfall(pC, mx=16, drng=7)
f.colorbar(im1, ax=ax1)
f.colorbar(im2, ax=ax2)
#f.colorbar()

#plt.show()

import IPython; IPython.embed()
