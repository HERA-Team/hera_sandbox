#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import numpy as np, aipy, capo, matplotlib.pyplot as plt, sys, glob
import md5
import oqe, os
import cProfile
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

SEPS = [(0,103), (0,111), (0,95)]
#SEPS += [(2,105), (1,83)]
#SEPS += [(0,79), (0,78)]
#SEPS += [(0,70),(0,71)]
#SEPS += [(1,107),(0,51)]
#SEPS += [(3,105),(3,106)]
#CH0,NCHAN = 90, 31
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
    dataDIR = '/Users/yunfanz/DATA/DATA128/DATA/'
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
    
#k1a,k1b,k1c = [(s,'xx',(0,103)) for s in sets]
#k2a,k2b,k2c = [(s,'xx',(0,111)) for s in sets]
#k3a,k3b,k3c = [(s,'xx',(0, 95)) for s in sets]
#ks = [k1a,k1b,k1c,k2a,k2b,k2c,k3a,k3b,k3c]
#k1a,k1b = [(s,'xx',(0,103)) for s in sets]
#k2a,k2b = [(s,'xx',(0,111)) for s in sets]
#k3a,k3b = [(s,'xx',(0, 95)) for s in sets]
#ks = [k1a,k1b,k2a,k2b,k3a,k3b]
ks = [(s,'xx',bl) for bl in SEPS for s in sets]
#k1a, = [(s,'xx',(0,103)) for s in sets]
#k2a, = [(s,'xx',(0,111)) for s in sets]
#k3a, = [(s,'xx',(0, 95)) for s in sets]
#ks = [k1a,k2a,k3a]
NK = len(ks)

def set_C(norm=3e-6):
    ds.clear_cache()
    Cs,iCs = {},{}
    for k in ks:
        #Cs[k] = sum([capo.oqe.cov(ds.x[k][:,400:],ds.w[k][:,400:])+norm*np.identity(NCHAN) for ki in ks if ki != k])
        #Cs[k] = sum([capo.oqe.cov(ds.x[ki][:,400:],ds.w[ki][:,400:])+norm*np.identity(NCHAN) for ki in ks if ki[2] != k[2]])
        Cs[k] = sum([oqe.cov(ds.x[k][:,400:],ds.w[k][:,400:])+norm*np.identity(NCHAN) for ki in ks if ki[2] != k[2]])
        #w = np.where(ds.w[ki] > 0, 1, 0)
        #Cs[k] = sum([capo.oqe.cov(ds.x[ki][:,400:],w[:,400:])+norm*np.identity(NCHAN) for ki in ks if ki[2] != k[2]])
        #Cs[k] = sum([ds.C(k)+norm*np.identity(NCHAN) for ki in ks if ki != k])
        #Cs[k] = sum([ds.C(k)+norm*np.identity(NCHAN) for ki in data if ki[2] != k[2]])
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
            qC = ds.q_hat(k1,k2)
            FC = ds.get_F(k1,k2)
            MC,WC = ds.get_MW(FC, mode='F^-1/2')
            pC = ds.p_hat(MC,qC)
        return pC * scalar

data_g, wgt_g = {},{}
for k in data:
    lst_g,data_g[k],wgt_g[k] = oqe.lst_grid(lsts[k[0]],data[k])
################################
ds = oqe.DataSet(data_g, wgt_g)

set_C(3e-6)
#pI,pW,pC = get_p(ks[0],ks[1])

for cnt,k in enumerate(ks):
    plt.subplot(NK,1,cnt+1)
    capo.plot.waterfall(ds.x[k], drng=3)
    plt.title(k)
    plt.colorbar()
plt.savefig('timeseries.png')

for cnt,k in enumerate(ks):
    plt.subplot(NK,1,cnt+1)
    pC = get_p(k,k,'C')
    plt.title(k)
    capo.plot.waterfall(pC, mx=16, drng=7)
    plt.colorbar()
plt.savefig('pspc.png')

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

import IPython; IPython.embed()
