#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys

MX = -1
DRNG = 2.5
seps = [
    #'<0,103> xx', #  1
    #'<1,4> xx',   #  2
    #'<0,101> xx', #  3
    #'<0,62> xx',  #  4
    #'<0,100> xx', #  5
    #'<1,13> xx',  #  6
    #'<1,70> xx',  #  7
    #'<1,56> xx',  #  8
    #'<1,71> xx',  #  9
    #'<1,59> xx',  # 10
    '<0,97> xx',  # 11
    '<12,43> xx', # 12
    '<9,71> xx',  # 13
    '<9,59> xx',  # 14
    '<57,64> xx', # 15
]

conj = ['<0,97> xx', '<12,43> xx', '<57,64> xx']

d,w = {}, {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    for sep in seps:
        if sep in conj: d[sep] = d.get(sep,[]) + [npz[sep].conj()]
        else: d[sep] = d.get(sep,[]) + [npz[sep]]
d_cat = n.concatenate([n.concatenate(d[k]) for k in seps], axis=1)
Cov = n.dot(d_cat.T, d_cat.conj())
p.subplot(121); C.arp.waterfall(Cov, drng=3); p.colorbar()
p.subplot(122); C.arp.waterfall(Cov, mode='phs'); p.colorbar(); p.show()
_d,_w = {},{}
for i,k in enumerate(d.keys()):
    d[k] = n.concatenate(d[k])
    fqs = n.linspace(.1,.2,d[k].shape[-1])
    #d = d[:,90:129]
    d[k],fqs = d[k][:,14:180], fqs[14:180]
    w[k] = n.where(d[k] == 0, 0, 1.)
    #w = n.where(n.isnan(d), 0, w)
    #w = n.where(w < 1e-8, 0, w)
    #d = n.where(w > 0, d, 0)
    p.figure(1); p.subplot(2,len(seps),i+1); C.arp.waterfall(d[k], mx=0, drng=4)
    window = a.dsp.gen_window(fqs.size, 'blackman-harris');
    window.shape = (1,window.size)
    _d[k] = n.fft.ifft(d[k]*window)
    _w[k] = n.fft.ifft(w[k]*window)
    p.figure(2)
    p.subplot(len(seps), 3, 3*i+1); C.arp.waterfall(_d[k], drng=4); p.colorbar()
    p.subplot(len(seps), 3, 3*i+2); C.arp.waterfall(_w[k], drng=4)
    area = n.ones(d[k].shape[1], dtype=n.int)
    #area[25:175] = 0
    #area[13:190] = 0
    #area[30:-31] = 0
    area[22:-21] = 0
    #area[9:-8] = 0
    dc = d[k].copy()
    _d2 = []
    for j in range(d[k].shape[0]):
        g = n.sqrt(n.average(w[k][j]**2))
        if g == 0: continue
        _dcl,info = a.deconv.clean(_d[k][j], _w[k][j], tol=1e-9, area=area, stop_if_div=False, maxiter=100)
        print k, j, info['term']
        dmdl = n.fft.fft(_dcl)
        _dcl += info['res'] / g
        _d2.append(_dcl)
        dc[j] -= dmdl * w[k][j]
    _d[k] = n.array(_d2)
    d[k] = dc
    print _d[k].shape
    p.figure(2); p.subplot(len(seps), 3, 3*i+3); C.arp.waterfall(_d[k], drng=4); p.colorbar()
    p.figure(1); p.subplot(2,len(seps),len(seps)+i+1); C.arp.waterfall(dc, mx=0, drng=4)
p.show()


d = n.concatenate([d[k] for k in seps], axis=1)
C.arp.waterfall(n.dot(d.T,d.conj()), drng=6)
p.show()

Cov = n.dot(d.T, d.conj())
#Cov = n.dot(d.T, d.conj())
print Cov.shape
U,S,V = n.linalg.svd(Cov)
p.semilogy(S); p.show()
#iS = n.where(S > 1e-6, 1./S, 0.)
iS = 1. / (S+1e-6)
_C = n.einsum('ij,j,jk', V.T, iS, U.T)
p.subplot(121); C.arp.waterfall(Cov, drng=5); p.colorbar(shrink=.5)
p.subplot(122); C.arp.waterfall(_C, drng=5); p.colorbar(shrink=.5); p.show()

_Cd = n.dot(_C, d.T).T
#_Cdc = n.dot(_C, d.T).T

p.subplot(221); C.arp.waterfall(d, mode='log', mx=MX, drng=DRNG); p.colorbar(shrink=.5)
p.subplot(222); C.arp.waterfall(d, mode='phs'); p.colorbar(shrink=.5)
p.subplot(223); C.arp.waterfall(_Cd, drng=DRNG, mode='log'); p.colorbar(shrink=.5)
p.subplot(224); C.arp.waterfall(_Cd, mode='phs'); p.colorbar(shrink=.5)
p.show()


import IPython; IPython.embed()
