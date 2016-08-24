#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys

MX = -1
DRNG = 2.5
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
    '<0,97> xx', # 11
    '<12,43> xx', # 12 
    '<57,64> xx' # 15
]

def dot(a1, a2, med=False):
    rv = n.zeros((a1.shape[0],a2.shape[1]), dtype=a1.dtype)
    for i in xrange(rv.shape[0]):
        for j in xrange(rv.shape[1]):
            a1a2 = a1[i] * a2[:,j]
            valid = n.where(n.abs(a1a2) > 0, 1, 0).sum()
            if valid == 0: valid = 1
            if med: rv[i,j] = n.median(a1a2)
            else: rv[i,j] = n.sum(a1a2) / valid
    return rv

d,w = {}, {}
chisq = []
for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    for sep in seps:
        if sep in conj: d[sep] = d.get(sep,[]) + [npz[sep].conj()]
        else: d[sep] = d.get(sep,[]) + [npz[sep]]
    csq = npz['chisq xx']
    U,S,V = n.linalg.svd(csq, full_matrices=False); S[3:] *= 0
    csq_ = n.dot(U, n.dot(n.diag(S), V))
    chisq.append(n.where(csq > 0, csq - csq_, 0))
    #chisq.append(n.where(csq > 0, csq - n.median(csq, axis=0), 0))
for sep in seps: d[sep] = n.concatenate(d[sep])
chisq = n.concatenate(chisq)
fq1 = n.linspace(.1,.2,d.values()[0].shape[1])
ch1 = n.arange(fq1.size)
sdf = fq1[1]-fq1[0]

if False: # full cross-mult
    d_cat = n.concatenate([d[k] for k in seps], axis=1)
    print d_cat.shape
    Cov = n.dot(d_cat.T, d_cat.conj())
    mCov = mdot(d_cat.T, d_cat.conj())
    p.subplot(121); C.arp.waterfall(Cov, mx=0, drng=3); p.colorbar()
    #p.subplot(122); C.arp.waterfall(Cov, mode='phs'); p.colorbar(); p.show()
    p.subplot(122); C.arp.waterfall(mCov, mx=0, drng=3); p.colorbar(); p.show()

d1 = d[seps[1]]
Cov1 = dot(d1.T, d1.conj())
#Cov2 = dot(d1.T, d1.conj(), med=True)

#p.subplot(121); C.arp.waterfall(Cov1, mx=0, drng=3); p.colorbar()
#p.subplot(122); C.arp.waterfall(Cov2, mx=0, drng=3); p.colorbar(); p.show()

U1,S1,V1 = n.linalg.svd(Cov1)
#U2,S2,V2 = n.linalg.svd(Cov2)
iCov1 = n.dot(V1.T.conj(), n.dot(n.diag(1/(S1+1e-6)), U1.T.conj()))
#iCov2 = n.dot(V2.T.conj(), n.dot(n.diag(1/(S2+.01)), U2.T.conj()))

ad1 = n.abs(d1)
d1C1 = n.where(ad1 > 0, n.dot(iCov1,d1.T).T, 0)
#d1C2 = n.where(ad1 > 0, n.dot(iCov2,d1.T).T, 0)

p.subplot(131); C.arp.waterfall(d1, drng=3); p.colorbar()
p.subplot(132); C.arp.waterfall(d1C1, drng=3); p.colorbar()
#p.subplot(133); C.arp.waterfall(d1C2, drng=3); p.colorbar()
p.subplot(133); C.arp.waterfall(chisq, drng=5); p.colorbar()
p.show()

import IPython; IPython.embed()

