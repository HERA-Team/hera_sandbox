#! /usr/bin/env python

import numpy as np, aipy, capo, pylab as plt, sys, glob

def dB(sig): return 10*np.log10(np.abs(np.average(sig.real, axis=1)))

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

SEP = SEPS[1]
#CH0,NCHAN = 90, 31
CH0,NCHAN = 110, 51

meta, gains, vismdl, xtalk = capo.omni.from_npz(sys.argv[1:], verbose=True)
data,wgts = {}, {}
for pol in vismdl:
    #for bl in vismdl[pol]:
    for bl in SEPS:
        k = ('day0',pol,bl)
        data[k] = vismdl[pol][bl][:,CH0:CH0+NCHAN]
        if bl in CONJ: data[k] = data[k].conj()
        wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1)
ds = capo.oqe.DataSet(data, wgts)
k1 = ('day0','xx',(0,103))

'''
Cs,iCs = {},{}
for k in data:
    #Cs[k] = ds.C(k)
    #Cs[k] = sum([ds.C(ki)+0*np.identity(NCHAN) for ki in dat])
    Cs[k] = sum([ds.C(ki)+3e-6*np.identity(NCHAN) for ki in dat if ki != k])
    #Cs[k] = sum([ds.C(ki)+1e-6*np.identity(NCHAN) for ki in dat if ki != k])
    #Cs[k] = sum([ds.C(ki)+1e-4*np.identity(NCHAN) for ki in dat if ki != k])
    #Cs[k] = sum([ds.C(ki)+0*np.identity(NCHAN) for ki in dat if ki != k])
    #ds.set_C({k:Cs[k]+1e-1*np.identity(NCHAN)}) # regularize a bit with some diagonal
    ds.set_C({k:Cs[k]})
    iCs[k] = ds.iC(k)
'''

#ds.set_data(dat_cut)
#ds.set_data(dat)
#tau = np.fft.fftshift(dly)
win1 = aipy.dsp.gen_window(NCHAN, 'blackman-harris'); win1.shape = (-1,1)
win2 = aipy.dsp.gen_window(NCHAN, 'blackman-harris')**1.5; win2.shape = (-1,1)
#for ki in iCs:
for ki in data:
    print ki
    qI = ds.q_hat(ki,ki,use_cov=False)
    FI = ds.get_F(ki,ki,use_cov=False)
    MI,WI = ds.get_MW(FI, mode='I')
    pI = ds.p_hat(MI,qI)
    pW1 = 1.6*2*np.abs(np.fft.fftshift(np.fft.ifft(win1*data[ki].T, axis=0), axes=0))**2
    pW2 = 2.4*2*np.abs(np.fft.fftshift(np.fft.ifft(win2*data[ki].T, axis=0), axes=0))**2
    #plt.figure(1)
    #plt.plot(tau, dB(pI), 'b', label='I')
    #plt.plot(tau, dB(pW1), 'g', label='W')
    #plt.plot(tau, dB(pW2), 'g', label='W')
    #plt.plot(tau, dB(pI_eor), 'k', label='E')
    #ds.set_iC({ki:iCs[ki]})
    qC = ds.q_hat(ki,ki)
    FC = ds.get_F(ki,ki)
    MC,WC = ds.get_MW(FC, mode='F^-1/2')
    pC = ds.p_hat(MC,qC)
    for cnt,pk in enumerate([pI,pW1,pW2,pC]):
        #plt.figure(1)
        plt.subplot(5,1,cnt+1); capo.plot.waterfall(pk, mx=-2, drng=6), plt.colorbar()
        #plt.figure(2); plt.plot(dB(pk))
    plt.subplot(5,1,5); capo.plot.waterfall(ds.x[ki], drng=3), plt.colorbar()
    plt.show()

import IPython; IPython.embed()
