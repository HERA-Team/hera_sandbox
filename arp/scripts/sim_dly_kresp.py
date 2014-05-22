#! /usr/bin/env python

import numpy as n, aipy as a, pylab as p, capo as C

F_c = .150
B = .02

#urange = (535,550)
#urange = (541,544)
urange = (542,543)

b = n.array([100., 0., 0]) # topocentric ns
L = n.linspace(-1,1,1024)
_L = n.fft.fftfreq(L.size, L[1]-L[0])
_L = n.fft.fftshift(_L)
fq = n.arange(F_c-B/2, F_c+B/2, .1/203)[:-1]
_fq = n.fft.fftfreq(fq.size, fq[1]-fq[0])
_fq = n.fft.fftshift(_fq)
aa = a.cal.get_aa('psa898_v003', fq)
z = n.sqrt(1-L**2)
bm = n.where(z > 0, aa[0].bm_response((L,n.zeros_like(L), z)), 0)

w_fq = a.dsp.gen_window(fq.size, 'blackman-harris') # doesnt work?
w_fq.shape = (w_fq.size,1)
bm *= w_fq
#bm = n.ones_like(bm)

i,j = n.indices((fq.size,L.size))
Ls,fs = L[j], fq[i]
us = b[0] * fs

z2x_1 = []
z2x_2 = []
for tau in _fq:
    d = bm * n.exp(-2j*n.pi*(_L[urange[0]:urange[1]][0]*Ls + tau*fs))
    d = n.fft.fftshift(d)
    d /= (d[0,0] / n.abs(d[0,0])) # set the absolute phase, just for consistency
    _d = n.fft.ifft2(d)
    _d = n.fft.fftshift(_d)
    z2x_2.append(_d[:-1:2,urange[0]:urange[1]].flatten()) # WORKS
    z2x_1.append(_d[:,urange[0]:urange[1]].flatten()) # DOESN'T WORK

    #p.subplot(121)
    ##C.arp.waterfall(d, extent=(L[0],L[-1],fq[0],fq[-1]), mode='real')
    #C.arp.waterfall(d, mode='real')
    #p.colorbar(shrink=.5)
    #p.xlabel('L')
    #p.ylabel(r'$\nu$')
    #p.subplot(122)
    ##C.arp.waterfall(_d, extent=(_L[0],_L[-1],_fq[0],_fq[-1]), drng=10)
    ##C.arp.waterfall(_d, drng=10)
    #C.arp.waterfall(_d, mode='real')
    #p.colorbar(shrink=.5)
    #p.xlabel('u')
    #p.ylabel(r'$\eta$')
    #p.show()

z2x = n.array(z2x_1); print 'z2x:', z2x.shape, z2x.dtype
#z2x = n.array(z2x_2); print 'z2x:', z2x.shape, z2x.dtype
#z2x = n.identity(z2x.shape[0])
SH = (z2x.shape[1],urange[1]-urange[0])
#for i in range(z2x.shape[0]): p.plot(z2x[i].real)
#p.show()


def random_phase(shape):
    return n.exp(2j*n.pi*n.random.uniform(size=shape))

z = random_phase(z2x.shape[-1])
amp = n.sin(n.arange(z.size).astype(n.float64)*2*n.pi/z.size)
#z *= amp
z = amp
x = n.einsum('ij,j', z2x, z) # worried about a conjugation here in z2x: doesn't matter b/c it gets squared
p_true = n.abs(z)**2
NOISE_AMP = 0.0
x += random_phase(x.shape) * NOISE_AMP * n.random.normal(size=x.shape)
print 'x:', x.shape, x.dtype

def dagger(M):
    axes = range(M.ndim)
    axes[-2],axes[-1] = axes[-1],axes[-2]
    return n.transpose(n.conj(M), axes)

#S = z2x.dot(z2x.T.conj()); print 'S:', S.shape, S.dtype
S = n.einsum('ij,kj', z2x, z2x.conj()); print 'S:', S.shape, S.dtype
N = n.identity(S.shape[0], dtype=n.complex128) * NOISE_AMP**2
Cov = S + N; print 'C:', Cov.shape, Cov.dtype
u,s,v = n.linalg.svd(Cov)
# XXX the number of eigenmodes preserved here in proj and Cinv affects normalization later
print n.around(n.abs(s), 4)
#Cinv = dagger(v).dot(n.diag(n.where(s>1e-15,1./s,0))).dot(dagger(u))
Cinv = n.einsum('jh,ji,ki', v.conj(), n.diag(n.where(s>1e-15,1./s,0)), u.conj())
p.subplot(131); C.arp.waterfall(Cov, drng=3); p.colorbar(shrink=.5)
p.subplot(132); C.arp.waterfall(Cinv, drng=3); p.colorbar(shrink=.5)
p.subplot(133); C.arp.waterfall(n.dot(Cov,Cinv), drng=3); p.colorbar(shrink=.5)
p.show()
z2x_t = n.transpose(z2x); z2x_t.shape += (1,); print 'z2x_t:', z2x_t.shape, z2x_t.dtype
Qa = z2x_t * dagger(z2x_t); print 'Qa:', Qa.shape, Qa.dtype
Ea = 0.5 * n.einsum('ij,ajk,kl', Cinv, Qa, Cinv); print 'Ea:', Ea.shape, Ea.dtype
EaC = n.einsum('aij,jk', Ea, Cov); print 'EaC:', EaC.shape, EaC.dtype
Sab = 2*n.einsum('aij,bji', EaC, EaC); print 'Sab:', Sab.shape, Sab.dtype
print n.diag(Sab)
Wab = n.einsum('aij,bji', Ea, Qa); print 'Wab:', Wab.shape, Wab.dtype
p.subplot(121); C.arp.waterfall(Wab, drng=3); p.colorbar(shrink=.5)
p.subplot(122); C.arp.waterfall(Sab, drng=3); p.colorbar(shrink=.5)
p.show()
if False:# XXX this renormalization doesn't reproduce input signal for singular C
    norm = Wab.sum(axis=-1); norm.shape += (1,1); print 'norm:', norm.shape
    Ea /= norm
else: # this renormalization works, regardless of the signularity/projection out of modes in Cinv
    Ea /= 0.5
pa = n.einsum('j,ajk,k', x.conj(), Ea, x); print 'pa:', pa.shape, pa.dtype
ba = n.einsum('aji,ij', Ea, N); ba.shape = pa.shape; print 'ba:', ba.shape, ba.dtype
pa -= ba
pa, ba, p_true = pa.flatten(), ba.flatten(), p_true.flatten()
p.subplot(131); p.plot(n.abs(z.flatten())**2); p.title('z')
p.subplot(132); p.plot(n.abs(x.flatten())**2); p.title('x')
p.subplot(133)
p.title('pa')
p.plot(p_true.real)
p.plot(pa.real, '.')
#p.plot(pa.real+n.diag(Sab), '.')
#p.plot(pa.real-n.diag(Sab), '.')
p.plot(ba.real)
#p.plot(n.arange(pa.size-1) + .5, 0.5*(pa[1:] + pa[:-1]).real, '+')
p.show()
#QC = n.einsum('ij,ajk', Cinv, Qa); print 'QC:', QC.shape, QC.dtype
#F = 0.5 * n.einsum('aji,bij', QC, QC); print 'F:', F.shape, F.dtype
#Fu,Fs,Fv = n.linalg.svd(F)
#print n.around(n.abs(Fs), 4)
# The whole normalization problem with projecting out modes comes down to the change
# in amplitude of eigenvalues in the Fisher matrix
#if False: F = Fu.dot(0.5*n.diag(n.ones_like(Fs)).dot(Fv)); print F.shape, F.dtype
#Finv = dagger(Fv).dot(n.diag(n.where(Fs>1e-3,1./Fs,0))).dot(dagger(Fu))
#Finv_sqrt = dagger(Fv).dot(n.diag(n.where(Fs>1e-3,1./n.sqrt(Fs),0))).dot(dagger(Fu))
#M = n.identity(F.shape[0], dtype=n.complex128); print 'M:', M.shape, M.dtype
#M = Finv; print 'M:', M.shape, M.dtype
#M = Finv_sqrt; print 'M:', M.shape, M.dtype
#W = n.dot(M, F); print 'W:', W.shape, W.dtype
#norm = W.sum(axis=-1); norm.shape += (1,); print 'norm:', norm.shape
#M /= norm 
#W = n.dot(M, F)
#pa = n.einsum('ba,aij', M, qa-na); print 'pa:', pa.shape, pa.dtype




