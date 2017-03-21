#! /usr/bin/env python
import aipy as a, capo as C, numpy as n, pylab as p
import sys

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

times,dat,flg = C.arp.get_dict_of_uv_data(sys.argv[1:], '0_7', 'xx', verbose=True)
d = dat.values()[0].values()[0]
print 'Computing covariance'
dd = n.cov(d.T)
print 'SVD'
U,S,V = n.linalg.svd(dd.conj())
_S = n.ones_like(S)
#_S = S.copy()
#_S[:30] = 0; _S[30:] = 1
#_S = n.where(S > 1e20, 0, 1)
#V = n.where(n.abs(V) > 1e-3, V, 0)
#U = n.where(n.abs(U) > 1e-3, U, 0)
#_dd = n.dot(V.T, n.dot(n.diag(_S), U.T))
_dd = n.diag(_S)

d_ = n.dot(_dd, d.T).T

p.subplot(321)
C.arp.waterfall(d, mx=8, drng=4)
p.colorbar()

p.subplot(322)
C.arp.waterfall(d_, mx=8, drng=4)
p.colorbar()

p.subplot(323)
C.arp.waterfall(d, mode='phs')
p.colorbar()

p.subplot(324)
C.arp.waterfall(d_, mode='phs')
p.colorbar()

p.subplot(325)
C.arp.waterfall(dd)
#p.semilogy(S)

p.subplot(326)
C.arp.waterfall(V)
p.colorbar()

p.show()
