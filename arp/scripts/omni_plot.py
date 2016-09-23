#! /usr/bin/env python
import numpy as np, capo, aipy as ap
import pylab as plt
import sys, glob

seps = [
    '<0,103> xx', #  1
    '<1,4> xx',   #  2
    '<0,101> xx', #  3
    '<0,62> xx',  #  4
    '<0,100> xx', #  5
    '<1,13> xx',  #  6
    '<1,70> xx',  #  7
    '<1,56> xx',  #  8
    '<1,71> xx',  #  9
    '<1,59> xx',  # 10
    '<0,97> xx',  # 11
    '<12,43> xx', # 12
    '<9,71> xx',  # 13
    '<9,59> xx',  # 14
    '<57,64> xx', # 15
]

conj = [
    '<0,103> xx', #  1
    '<1,4> xx',   #  2
    '<0,101> xx', #  3
    '<0,62> xx',  #  4
    '<0,100> xx', #  5
    '<0,97> xx',  # 11
    '<12,43> xx', # 12
    '<57,64> xx'] # 15

SEP = seps[1]

def cov(d1, w1, d2=None, w2=None):
    if d2 is None: d2,w2 = d1.conj(),w1
    d1sum,d1wgt = d1.sum(axis=1), w1.sum(axis=1)
    d2sum,d2wgt = d2.sum(axis=1), w2.sum(axis=1)
    x1,x2 = d1sum / np.where(d1wgt > 0,d1wgt,1), d2sum / np.where(d2wgt > 0,d2wgt,1)
    x1.shape = (-1,1)
    x2.shape = (-1,1)
    d1x = d1 - x1
    d2x = d2 - x2
    C = np.dot(d1x,d2x.T)
    W = np.dot(w1,w2.T)
    return C / np.where(W > 1, W-1, 1)

def get_data(filelist, seps, conj, verbose=True):
    d,w = {}, {}
    lsts, jds = [], []
    for filename in filelist:
        if verbose: print 'Reading', filename
        npz = np.load(filename)
        for sep in seps: d[sep] = d.get(sep,[]) + [npz[sep]]
        lsts.append(npz['lsts'])
        jds.append(npz['jds'])
    for sep in seps:
        d[sep] = np.concatenate(d[sep])
        w[sep] = np.where(d[sep] != 0, 1., 0)
    for sep in conj:
        try: d[sep] = d[sep].conj()
        except(KeyError): pass
    return np.concatenate(lsts), np.concatenate(jds), d, w

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

for SEP in seps[:1]:
    binsize = 0.00230769985489
    lst,d,w = {},{},{}
    #for jd in ['980','981','982','983','985','986','987']:
    for jd in ['980','981','982']:
        files = glob.glob('*%s.*.xx.uvcRRE.npz' % jd)
        lst_jd,_,d_jd,w_jd = get_data(files, [SEP], conj)
        lst[jd],d[jd],w[jd] = rebin_lst(binsize, lst_jd, d_jd[SEP], w_jd[SEP])

    # Align lst waterfalls for each day (XXX may want to do this pairwise)
    lstmin = max([lst[jd][0] for jd in lst.keys()])
    lstmax = min([lst[jd][-1] for jd in lst.keys()])

    for jd in lst.keys():
        print jd, lst[jd][0],lst[jd][-1], lstmin, lstmax
        i = np.argwhere(lst[jd] == lstmin)[0]
        j = np.argwhere(lst[jd] == lstmax)[0]
        d[jd], w[jd], lst[jd] = d[jd][i:j], w[jd][i:j], lst[jd][i:j]

    dsum = sum([d[jd]*w[jd] for jd in lst.keys()])
    dwgt = sum([w[jd] for jd in lst.keys()])
    davg = dsum / np.where(dwgt > 0, dwgt, 1)
    dmed = np.median([d[jd] for jd in lst.keys()], axis=0)

    if True:
        MX,DRNG = 0,3
        for i, jd in enumerate(lst.keys()):
            plt.subplot(1,len(lst)+1,i+1)
            mask = np.where(w[jd] > 0, 1, 0)
            #ddi = d[jd] - davg*mask
            ddi = d[jd]
            C = cov(ddi.T,w[jd].T)
            U,S,V = np.linalg.svd(C)
            cutoff = np.median(S)
            _C = np.dot(V.T.conj(), np.dot(np.diag(1/(S+cutoff)), U.T.conj()))
            #capo.arp.waterfall(d[jd], mx=0, drng=3)
            capo.arp.waterfall(np.dot(_C,d[jd].T).T, drng=DRNG)
            #capo.arp.waterfall(d[jd]-dmed*mask, mx=0, drng=3)
            #capo.arp.waterfall(d[jd]-dmed*mask, mode='phs')
        plt.subplot(1,len(lst)+1,len(lst)+1)
        capo.arp.waterfall(davg, mx=MX, drng=DRNG)
        #capo.arp.waterfall(dmed, mx=0, drng=3)
        #capo.arp.waterfall(dmed, mode='phs')
        plt.show()

        for i, jdi in enumerate(lst.keys()):
            for j,jdj in enumerate(lst.keys()):
                plt.subplot(len(lst),len(lst),i*len(lst)+j+1)
                capo.arp.waterfall(cov(d[jdi].T,w[jdi].T,d[jdj].T.conj(),w[jdj].T), mx=0, drng=4)
        plt.show()

        for i, jdi in enumerate(lst.keys()):
            mask = np.where(w[jdi] > 0, 1, 0)
            di = d[jdi] - dmed*mask
            for j,jdj in enumerate(lst.keys()):
                mask = np.where(w[jdj] > 0, 1, 0)
                dj = d[jdj] - dmed*mask
                plt.subplot(len(lst),len(lst),i*len(lst)+j+1)
                #capo.arp.waterfall(np.dot(di.T, dj.conj()), drng=3)
                capo.arp.waterfall(cov(di.T,w[jdi].T,dj.T.conj(),w[jdj].T), mx=-4, drng=4)
        plt.show()

    dm = {}
    for jd in d:
        dm[jd] = d[jd] * dmed.conj()
    mm = dmed * dmed.conj()
    mm = mm[:850]
    #mmi = mm

    for i,jd in enumerate(dm):
        #plt.subplot(1,4,i+1)
        print dm[jd].shape
        #dmi = dm[jd]
        dmi = dm[jd][:850]
        dmi = dmi.reshape((dmi.shape[0]/50, -1, dmi.shape[1]))
        mmi = mm.reshape((mm.shape[0]/50, -1, mm.shape[1]))
        #capo.arp.waterfall(np.sum(dmi,axis=1)/np.sum(mmi,axis=1), mx=1.1, drng=.2, mode='real')
        dmi_mmi = np.sum(dmi,axis=1)/np.sum(mmi,axis=1)
        plt.plot(np.median(dmi_mmi,axis=1))
plt.show()

import IPython; IPython.embed()
