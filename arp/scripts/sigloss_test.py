#! /usr/bin/env python
import numpy as np, pylab as plt
import capo, aipy
import capo.oqe as oqe
import sys

CH0,NCHAN = 30, 61
NSAMP = 120

for i in xrange(10):
    e = oqe.noise(size=(NCHAN,NSAMP))
    v = oqe.noise(size=(NCHAN,NSAMP))
    r = e + v
    
    k = ('even',(0,1),'I')
    k1 = ('e',(0,1),'I')
    k2 = ('v',(0,1),'I')
    ds = oqe.DataSet(dsets={k:r.T[:-20]})
    print i, np.linalg.cond(ds.C(k))
    iC_r = ds.iC(k)
    ds.set_data({k:e.T}); q_e = ds.q_hat(k,k)
    ds.set_data({k:v.T}); q_v = ds.q_hat(k,k)
    ds.set_data({k:r.T}); q_r = ds.q_hat(k,k)
    F = ds.get_F(k,k)
    ds.set_data({k1:e.T, k2:v.T})
    ds.set_iC({k1:iC_r, k2:iC_r})
    q_ev = ds.q_hat(k1,k2)
    (M,W) = ds.get_MW(F)
    p_e = ds.p_hat(M,q_e)
    p_v = ds.p_hat(M,q_v)
    p_r = ds.p_hat(M,q_r)
    p_ev = ds.p_hat(M,q_ev)
    print i
    print 'x', np.var(e), np.var(v), np.var(r)
    print 'q', np.mean(q_e), np.mean(q_v), np.mean(q_r), np.mean(q_ev)
    print 'p', np.mean(p_e), np.mean(p_v), np.mean(p_r), np.mean(p_ev)

plt.subplot(311); capo.plot.waterfall(e, mode='real'); plt.colorbar()
plt.subplot(312); capo.plot.waterfall(v, mode='real'); plt.colorbar()
plt.subplot(313); capo.plot.waterfall(r, mode='real'); plt.colorbar()
plt.show()


plt.subplot(411); capo.plot.waterfall(q_e , mode='real'); plt.colorbar()
plt.subplot(412); capo.plot.waterfall(q_v , mode='real'); plt.colorbar()
plt.subplot(413); capo.plot.waterfall(q_r , mode='real'); plt.colorbar()
plt.subplot(414); capo.plot.waterfall(q_ev, mode='real'); plt.colorbar()
plt.show()


plt.subplot(411); capo.plot.waterfall(p_e , mode='real'); plt.colorbar()
plt.subplot(412); capo.plot.waterfall(p_v , mode='real'); plt.colorbar()
plt.subplot(413); capo.plot.waterfall(p_r , mode='real'); plt.colorbar()
plt.subplot(414); capo.plot.waterfall(p_ev, mode='real'); plt.colorbar()
plt.show()

