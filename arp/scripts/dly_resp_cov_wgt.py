#! /usr/bin/env python
import numpy as np, pylab as plt
import capo, aipy
import capo.oqe as oqe

NCHAN = 100
Q = {}
for mode in xrange(NCHAN):
    Q[mode] = oqe.get_Q(mode, NCHAN)

#capo.plot.waterfall(Q[53], mode='real'); plt.show()

fg_ch = np.sin(np.linspace(0,2*np.pi,NCHAN)); fg_ch.shape = (-1,1)
fg_t = np.sin(np.linspace(0,2*np.pi, 1000)); fg_t.shape = (1,-1)
fg = fg_ch * fg_t

bp = aipy.dsp.gen_window(NCHAN, 'hamming'); bp.shape = (-1,1)

x = np.random.normal(size=(NCHAN,1000)) + fg
x *= bp
wgt = np.ones_like(x)

C = oqe.cov(x, wgt)
#U,S,V = np.linalg.svd(C)
#import IPython; IPython.embed()
_C = np.linalg.inv(C)
_Cx = np.dot(_C, x)
qC = np.abs(np.fft.ifft(_Cx, axis=0)); qC = np.average(qC, axis=1)
qI = np.abs(np.fft.ifft(x, axis=0)); qI = np.average(qI, axis=1)
print _Cx.shape
#q = {}
#for mode in xrange(NCHAN):
#    print mode
#    q[mode] = _Cx.conj() * np.dot(Q[mode], _Cx)
#    print q[mode].shape

#q = np.array([q[mode] for mode in xrange(NCHAN)])
#plt.subplot(211); capo.plot.waterfall(np.fft.fftshift(qI, axes=[0]), drng=3)
#plt.subplot(212); capo.plot.waterfall(np.fft.fftshift(qC, axes=[0]), drng=3)
plt.plot(np.fft.fftshift(qI))
plt.plot(np.fft.fftshift(qC))
plt.show()

#plt.figure(1)
#plt.subplot(211); capo.plot.waterfall(x, mode='real')
#plt.subplot(212); capo.plot.waterfall(_Cx, mode='real')
#plt.figure(2)
##plt.subplot(121); capo.plot.waterfall(C, drng=3)
#plt.subplot(121); capo.plot.waterfall(C, mode='real')
#plt.subplot(122); capo.plot.waterfall(_C, drng=3)
#plt.show()

