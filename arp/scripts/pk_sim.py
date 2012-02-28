#! /usr/bin/env python

import numpy as n, pylab as p

def k3pk(k): return n.where(k < .1, k**2*1e4, 1e2)
def dk3pk(k,dkx,dky,dkz): return k3pk(k) / (4*n.pi*k**3) * dkx * dky * dkz

#KMAX=2.
KMAX=.2
SIZE=128
TSYS=20
kx = n.arange(.1/SIZE, KMAX+2*KMAX/SIZE, 2*KMAX/SIZE)
ky = n.arange(.1/SIZE, KMAX+2*KMAX/SIZE, 2*KMAX/SIZE)
kz = n.arange(.1/SIZE, KMAX+2*KMAX/SIZE, 2*KMAX/SIZE)
ks = n.sqrt(kx.reshape((kx.size,1,1))**2 + ky.reshape((1,ky.size,1))**2 + kz.reshape((1,1,kz.size))**2)
_box = n.sqrt(dk3pk(ks, kx[1]-kx[0], ky[1]-ky[0], kz[1]-kz[0]))
print _box.shape
_box = n.concatenate([_box, _box[-2:0:-1]], axis=0)
ks = n.concatenate([ks, ks[-2:0:-1]], axis=0)
print _box.shape
_box = n.concatenate([_box, _box[:,-2:0:-1]], axis=1)
ks = n.concatenate([ks, ks[:,-2:0:-1]], axis=1)
print _box.shape
print _box[0,0]
_box *= n.exp(1j*n.random.uniform(0, 2*n.pi, size=_box.shape).astype(n.complex))
print _box[0,0]
_box = n.concatenate([_box, n.conj(_box)[::-1,::-1,-2:0:-1]], axis=2)
ks = n.concatenate([ks, ks[::-1,::-1,-2:0:-1]], axis=2)
_box[0,0,0] = n.abs(_box[0,0,0])
print _box[0,0]
print _box.shape
box = n.fft.fftn(_box)
box2 = box + n.random.normal(scale=TSYS, size=box.shape)

#T_eor = 0.010
#ks = 2*n.pi*n.arange(-.5, .5, 1./L)
#ks = n.concatenate([ks[ks.size/2:], ks[:ks.size/2]])
#box = n.random.normal(scale=T_eor, size=(L,L,L))
#_box = n.fft.ifftn(box)
print 'Parseval check:', n.average(n.abs(box)**2), n.sum(n.abs(_box)**2)

avg_over_pr = n.average(n.average(box, axis=0), axis=0)
_avg_over_pr = n.fft.ifft(avg_over_pr)
avg_over_pr2 = n.average(n.average(box2, axis=0), axis=0)
_avg_over_pr2 = n.fft.ifft(avg_over_pr2)
#k3pk_avg_over_pr = n.abs(_avg_over_pr)**2 / (2*n.pi/L)**3 * 4*n.pi*n.abs(ks)**3

p.subplot(211)
p.plot(box[box.shape[0]/2,box.shape[1]/2,:])
p.plot(avg_over_pr)
p.plot(box2[box2.shape[0]/2,box2.shape[1]/2,:])
p.plot(avg_over_pr2)

p.subplot(212)
p.loglog(ks[0,0], n.abs(_avg_over_pr)**2)
p.loglog(ks[0,0], n.abs(_avg_over_pr2)**2)

p.show()
