#! /usr/bin/env python

import numpy as n, aipy as a, pylab as p, capo as C

F_c = .150
B = .02
drng = 12

SCALE = 4.

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
w_fq *= n.sqrt(SCALE)
#w_fq /= w_fq.sum()
#w_fq = a.dsp.gen_window(fq.size, 'gaussian0.4') # works!
#w_fq = a.dsp.gen_window(fq.size, 'hamming') # works!
#w_fq = a.dsp.gen_window(fq.size, 'hanning') # doesn't work?
#w_fq = a.dsp.gen_window(fq.size, 'kaiser3') # almost works
#w_fq = n.ones_like(w_fq)
w_fq.shape = (w_fq.size,1)
bm *= w_fq
#bm = n.ones_like(bm)

i,j = n.indices((fq.size,L.size))
Ls,fs = L[j], fq[i]
us = b[0] * fs

z2x_1 = []
z2x_2 = []
for tau in _fq:
#for tau in _fq[::2]:
#for tau in n.fft.fftfreq(int(fq.size/2.5), fq[1]-fq[0]):
#for tau in n.fft.fftfreq(fq.size*8, fq[1]-fq[0]):
    #print tau
    #d = bm * n.exp(-2j*n.pi*(us*Ls + tau*fs))
    d = bm * n.exp(-2j*n.pi*(_L[urange[0]:urange[1]][0]*Ls + tau*fs))
    #dzero = n.zeros_like(d)
    #d = n.concatenate([dzero] * 5 + [d] + [dzero] * 5, axis=0)
    d = n.fft.fftshift(d)
    d /= (d[0,0] / n.abs(d[0,0])) # set the absolute phase, just for consistency
    _d = n.fft.ifft2(d)
    _d = n.fft.fftshift(_d)
    z2x_2.append(_d[:-1:2,urange[0]:urange[1]].flatten()) # WORKS
    #z2x_2.append(_d[:-1:4,urange[0]:urange[1]].flatten()) # WORKS
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
z2x_1 = n.array(z2x_1)
z2x_2 = n.array(z2x_2)
z2x = z2x_1; print 'z2x:', z2x.shape, z2x.dtype
#z2x = z2x_2; print 'z2x:', z2x.shape, z2x.dtype
SH = (z2x.shape[1],urange[1]-urange[0])
for i in range(z2x.shape[0]):
    #p.plot(n.abs(z2x[i]))
    p.plot(z2x[i].real)
    #print z2x[i].sum()
p.show()


def random_signal(shape):
    #amp = n.random.normal(size=shape)
    amp = 1.
    phs = n.exp(2j*n.pi*n.random.uniform(size=shape))
    return amp * phs

z = random_signal(z2x.shape[-1]); z.shape = (z.size,1)
#amp = n.zeros_like(z); amp[3:5,0] = 1.
#amp = n.zeros_like(z); amp[6:10,0] = 1.
#amp = n.zeros_like(z); amp[12:20:2,0] = 1.
amp = n.reshape(n.sin(n.arange(z.size).astype(n.float64)*2*n.pi/z.size),z.shape)
#z *= amp
z = amp
x = n.dot(z2x, z) # worried about a conjugation here in z2x: doesn't matter b/c it gets squared
NOISE_AMP = 0.
x += random_signal(x.shape) * NOISE_AMP * n.random.normal(size=x.shape)
#x = random_signal(x.shape)
print 'x:', x.shape, x.dtype

def dagger(M):
    axes = range(M.ndim)
    axes[-2],axes[-1] = axes[-1],axes[-2]
    return n.transpose(n.conj(M), axes)

S = z2x.dot(z2x.T.conj())
N = n.identity(S.shape[0], dtype=n.complex128) * NOISE_AMP**2
Cov = S + N; print 'C:', Cov.shape, Cov.dtype
u,s,v = n.linalg.svd(Cov)
# XXX the number of eigenmodes preserved here in proj and Cinv affects normalization later
print n.around(n.abs(s), 4)
Cinv = dagger(v).dot(n.diag(n.where(s>SCALE*1e-5,1./s,0))).dot(dagger(u))
print '|C|:', n.linalg.det(Cov)
#Cinv1 = n.linalg.inv(Cov)
#print 'C^-1:', Cinv.shape, Cinv.dtype, n.allclose(n.dot(Cov, Cinv), n.eye(Cov.shape[0]))
#p.subplot(121); C.arp.waterfall(Cinv, drng=3); p.colorbar(shrink=.5)
#p.subplot(122); C.arp.waterfall(Cinv1, drng=3); p.colorbar(shrink=.5)
#p.show()
p.subplot(131); C.arp.waterfall(Cov, drng=3); p.colorbar(shrink=.5)
p.subplot(132); C.arp.waterfall(Cinv, drng=3); p.colorbar(shrink=.5)
p.subplot(133); C.arp.waterfall(n.dot(Cov,Cinv), drng=3); p.colorbar(shrink=.5)
p.show()
z2x_t = n.transpose(z2x); z2x_t.shape += (1,); print 'z2x_t:', z2x_t.shape, z2x_t.dtype
Qa = z2x_t * dagger(z2x_t); print 'Qa:', Qa.shape, Qa.dtype
# Recompute Qa after projecting out modes in C that are singular
if False: # XXX is this necessary?  seems like not
    proj = u.dot(n.diag(n.where(s>SCALE*1e-15,1.,0))).dot(v) # project modes out of x that are singular
    p.subplot(121); C.arp.waterfall(Qa[0], drng=3)
    z2x_t = n.einsum('ji,aik',proj, z2x_t); print 'z2x_t:', z2x_t.shape, z2x_t.dtype
    Qa = z2x_t * dagger(z2x_t); print 'Qa:', Qa.shape, Qa.dtype
    p.subplot(122); C.arp.waterfall(Qa[0], drng=3)
    p.show()
xinv = n.dot(Cinv, x); print 'xinv:', xinv.shape, xinv.dtype
#p.plot(n.abs(x))
#p.plot(n.abs(z))
#p.plot(n.abs(xinv))
#p.show()
#qa = 0.5 * xinv.T.conj().dot(Qa.dot(xinv)); print 'qa:', qa.shape, qa.dtype
qa = 0.5 * n.einsum('ij,ajk', xinv.T.conj(), n.einsum('aji,ik', Qa, xinv)); print 'qa:', qa.shape, qa.dtype
# at this point, qa appears to be normalized correctly.  but pa is not
#qa = 0.5 * xinv.T.conj().dot(Qa.dot(xinv)); print 'qa:', qa.shape, qa.dtype
p.subplot(221); p.plot(n.abs(z.flatten())); p.title('z')
p.subplot(222); p.plot(n.abs(x.flatten())); p.title('x')
p.subplot(223); p.plot(n.abs(xinv.flatten())); p.title('x^-1')
p.subplot(224); p.plot(n.abs(qa.flatten())); p.title('qa')
# XXX Ea may need some work here: some things are off
Ea = 0.5 * n.einsum('ij,ajk', Cinv, n.einsum('aij,jk', Qa, Cinv)); print 'Ea:', Ea.shape, Ea.dtype
qa = n.einsum('ij,ajk', x.T.conj(), n.einsum('aij,jk', Ea, x)); print 'qa:', qa.shape, qa.dtype
p.subplot(224); p.plot(n.abs(qa.flatten())); p.title('qa')
p.show()
#Ea = Cinv.dot(Qa.dot(Cinv)).transpose([1,0,2]); print 'Ea:', Ea.shape, Ea.dtype
na = 0.5 * n.einsum('aji,ij', Ea, N); na.shape = qa.shape; print 'na:', na.shape, na.dtype
#QC = n.einsum('aij,jk', Qa, Cinv); print 'QC:', QC.shape, QC.dtype
QC = n.einsum('ij,ajk', Cinv, Qa); print 'QC:', QC.shape, QC.dtype
F = 0.5 * n.einsum('aji,bij', QC, QC); print 'F:', F.shape, F.dtype
Fu,Fs,Fv = n.linalg.svd(F)
print n.around(n.abs(Fs), 4)
# The whole normalization problem with projecting out modes comes down to the change
# in amplitude of eigenvalues in the Fisher matrix
if False: F = Fu.dot(0.5*n.diag(n.ones_like(Fs)).dot(Fv)); print F.shape, F.dtype
# F is not singular, so SCALE doesn't really matter here
Finv = dagger(Fv).dot(n.diag(n.where(Fs>SCALE*1e-3,1./Fs,0))).dot(dagger(Fu))
Finv_sqrt = dagger(Fv).dot(n.diag(n.where(Fs>SCALE*1e-3,1./n.sqrt(Fs),0))).dot(dagger(Fu))
M = n.identity(F.shape[0], dtype=n.complex128); print 'M:', M.shape, M.dtype
#M = Finv; print 'M:', M.shape, M.dtype
#M = Finv_sqrt; print 'M:', M.shape, M.dtype
W = n.dot(M, F); print 'W:', W.shape, W.dtype
norm = W.sum(axis=-1); norm.shape += (1,); print 'norm:', norm.shape
M /= norm 
W = n.dot(M, F)
pa = n.einsum('ba,aij', M, qa-na); print 'pa:', pa.shape, pa.dtype
W = n.einsum('aij,bji', Ea, Qa); print 'F:', F.shape, F.dtype
norm = W.sum(axis=-1); norm.shape += (1,1); print 'norm:', norm.shape
Ea /= norm
pa = n.einsum('ij,ajk', x.T.conj(), n.einsum('aij,jk', Ea, x)); print 'pa:', pa.shape, pa.dtype
#print W.sum(axis=-1) # verify normalization
p.subplot(121); C.arp.waterfall(F, drng=3); p.colorbar(shrink=.5)
#p.subplot(122); C.arp.waterfall(Finv, drng=3); p.colorbar(shrink=.5); p.show()
p.subplot(122); C.arp.waterfall(W, drng=3); p.colorbar(shrink=.5); p.show()
#p.subplot(221); C.arp.waterfall(n.reshape(W[0],SH) , drng=3)
#p.subplot(222); C.arp.waterfall(n.reshape(W[15],SH) , drng=3)
#p.subplot(223); C.arp.waterfall(n.reshape(W[30],SH) , drng=3)
#p.subplot(224); C.arp.waterfall(n.reshape(W[45],SH) , drng=3)
p.show()

#pa = M.dot(qa-na); print 'pa:', pa.shape, pa.dtype
#pa = M.dot(qa); print 'pa:', pa.shape, pa.dtype
pa.shape = SH
p_true = n.abs(z)**2
p_true.shape = pa.shape
p.plot(pa.flatten().real, '.')
p.plot(p_true.flatten().real)
#p.subplot(121); C.arp.waterfall(pa, mode='real'); p.colorbar(shrink=.5)
#p.subplot(122); C.arp.waterfall(p_true, mode='real'); p.colorbar(shrink=.5)
p.show()




