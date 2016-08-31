#! /usr/bin/env python

import numpy as np, aipy, capo, pylab as plt, sys, glob

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
    """ 
    rebin the lsts
    binsize: width of bins in seconds
    lsts: original lsts
    d: 
    """
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
    lstbin = (b0 + np.arange(bmax)) * binsize  #binned lsts
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
SEPS += [(2,105), (1,83)]
#SEPS += [(0,79), (0,78)]
#SEPS += [(0,70),(0,71)]
#SEPS += [(1,107),(0,51)]
#SEPS += [(3,105),(3,106)]
#CH0,NCHAN = 90, 31
CH0,NCHAN = 110, 51
dataDIR = '/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/'
sets = {
    #'day0' : sys.argv[1:],
    #'day0' : glob.glob('zen.2456714.*.xx.npz'),
    'day1' : glob.glob(dataDIR+'zen.2456715.52*.xx.npz'),
    'day2' : glob.glob(dataDIR+'zen.2456716.52*.xx.npz'),
}
data,wgts = {}, {}
lsts = {}
for s in sets:
    if not lsts.has_key(s):
        meta, gains, vismdl, xtalk = capo.omni.from_npz(sets[s], verbose=True)
        lsts[s] = meta['lsts']
    for pol in vismdl:
        #for bl in vismdl[pol]:
        for bl in SEPS:
            k = (s,pol,bl)
            data[k] = vismdl[pol][bl][:,CH0:CH0+NCHAN]
            if bl in CONJ: data[k] = data[k].conj()     #???
            wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1)              #set zero weight for where there's no data
#ds = capo.oqe.DataSet(data, wgts)
ind = {}
set1,set2 = sets.keys()[0], sets.keys()[-1]
lst_res = np.average(lsts[set1][1:] - lsts[set1][:-1])/2


################################################################
def lst_align(lsts1,lsts2,lstres,offset=0):
    i=0
    if lsts1[0]+offset>lsts2[0]:
        while lsts1[0]+offset>lsts2[i]: i += 1
        Npt = min(lsts1.size,lsts2.size-i)
        return np.arange(0,Npt),np.arange(i,Npt+i)
    else:
        while lsts1[i]+offset>lsts2[0]: i += 1
        Npt = min(lsts1.size-i,lsts2.size)
        return np.arange(i,Npt+i),n.arange(0,Npt)

ind[set1], ind[set2] = lst_align(lsts[set1], lsts[set2], lstres=lst_res)

for k in data.keys():
    (s,pol,bl) = k
    data[k] = data[k][ind[s]]     #aligning the lsts
    wgts[k] = wgts[k][ind[s]]
ds = capo.oqe.DataSet(data, wgts)
    
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
    ds.clear_cache()    #Cs are stored in cache to avoid recomputation
    Cs,iCs = {},{}
    for k in ks:
        #Cs[k] = sum([capo.oqe.cov(ds.x[k][:,400:],ds.w[k][:,400:])+norm*np.identity(NCHAN) for ki in ks if ki != k])

        import IPython; IPython.embed()
        Cs[k] = sum([capo.oqe.cov(ds.x[k][:,400:],ds.w[k][:,400:])+norm*np.identity(NCHAN) for ki in ks if ki[2] != k[2]])
        #Cs[k] = sum([ds.C(k)+norm*np.identity(NCHAN) for ki in ks if ki != k])
        #Cs[k] = sum([ds.C(k)+norm*np.identity(NCHAN) for ki in data if ki[2] != k[2]])
        ds.set_C({k:Cs[k]})
        iCs[k] = ds.iC(k)

#tau = np.fft.fftshift(dly)
win = aipy.dsp.gen_window(NCHAN, 'blackman-harris'); win.shape = (-1,1)   #???

def get_p(k1,k2,mode):
    assert(mode in 'IWC')   #IWC: identity, window, or covariance matrix
    if mode == 'I':
        qI = ds.q_hat(k1,k2,use_cov=False)
        FI = ds.get_F(k1,k2,use_cov=False)
        MI,WI = ds.get_MW(FI, mode='I')
        pI = ds.p_hat(MI,qI)
        return pI
    elif mode == 'W':
        pW = 1.6*2*np.fft.fftshift(np.fft.ifft(win*data[k1].T, axis=0) * np.fft.ifft(win*data[k2].T, axis=0).conj(), axes=0)
        return pW
    elif mode == 'C':
        qC = ds.q_hat(k1,k2)
        FC = ds.get_F(k1,k2)
        MC,WC = ds.get_MW(FC, mode='F^-1/2')
        pC = ds.p_hat(MC,qC)
        return pC

#set_C(1e-6)
set_C(0)  #??? norm=0
#pI,pW,pC = get_p(ks[0],ks[1])

for cnt,k in enumerate(ks):
    plt.subplot(NK,1,cnt+1)
    capo.plot.waterfall(ds.x[k], drng=3)
    plt.colorbar()
plt.show()

for cnt,k in enumerate(ks):
    plt.subplot(NK,1,cnt+1)
    pC = get_p(k,k,'C')
    plt.title(k[0])
    capo.plot.waterfall(pC, mx=-2, drng=7)
    plt.colorbar()
plt.show()

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

#import IPython; IPython.embed()
