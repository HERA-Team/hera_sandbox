#! /usr/bin/env python

import numpy as n, aipy as a, pylab as p, capo as C

ORTH = True
F_c = .150
B = .02
drng = 12

#b = n.array([100., 0., 0]) # topocentric ns
b = n.array([0., 0., 0]) # topocentric ns
L = n.linspace(-1,1,1024)
_L = n.fft.fftfreq(L.size, L[1]-L[0])
_L = n.fft.fftshift(_L)
fq = n.arange(F_c-B/2, F_c+B/2, .1/203)
_fq = n.fft.fftfreq(fq.size, fq[1]-fq[0])
_fq = n.fft.fftshift(_fq)
aa = a.cal.get_aa('psa898_v003', fq)
z = n.sqrt(1-L**2)
bm = n.where(z > 0, aa[0].bm_response((L,n.zeros_like(L), z)), 0)

w_fq = a.dsp.gen_window(fq.size, 'blackman-harris')
w_fq.shape = (w_fq.size,1)

i,j = n.indices((fq.size,L.size))
Ls,fs = L[j], fq[i]
us = b[0] * fs
tau = 100.

bm *= w_fq
bm = 1
d = bm * n.exp(-2j*n.pi*(us*Ls + tau*fs))
_d = n.fft.ifft2(d)
_d = n.fft.fftshift(_d)

if ORTH:
    k_u = -b[0] * F_c * C.pspec.dk_du(C.pspec.f2z(F_c)) # XXX shouldn't need minus sign here
    k_eta = -tau * C.pspec.dk_deta(C.pspec.f2z(F_c)) # XXX also shouldn't need minus sign here
    d_orth = n.zeros_like(d)
    r9 = 6595. # h^-1 Mpc to z=9. = 9421 * .7
    print k_u, k_eta
    def r(f): return r9 - C.pspec.dL_df(9.) * (f - C.pspec.z2f(9.))
    r150 = r(F_c)
    for i,f in enumerate(fq):
        r_i = r(f)
        x = r_i * L
        th = n.arcsin(L)
        y = r_i * n.cos(th)
        y -= r150
        p.plot(x,y)
        d_orth[i] = n.exp(-1j*(x*k_u + y*k_eta))
    d_orth = bm * d_orth
    _d_orth = n.fft.ifft2(d_orth)
    _d_orth = n.fft.fftshift(_d_orth)
    p.show()

if ORTH: p.subplot(221)
else: p.subplot(121)
p.imshow(n.real(d), extent=(L[0],L[-1],fq[0],fq[-1]), origin='lower', aspect='auto', interpolation='nearest')
p.colorbar(shrink=.5)
p.xlabel('L')
p.ylabel(r'$\nu$')

if ORTH: p.subplot(222)
else: p.subplot(122)
plt_d = n.log10(n.abs(_d)**2)
#mx = plt_d.max()
mx = 0
p.imshow(plt_d, extent=(_L[0],_L[-1],_fq[0],_fq[-1]), vmax=mx, vmin=mx-drng, origin='lower', aspect='auto', interpolation='nearest')
p.colorbar(shrink=.5)
p.xlabel('u')
p.ylabel(r'$\eta$')

if ORTH:
    p.subplot(223)
    p.imshow(n.real(d_orth), extent=(L[0],L[-1],fq[0],fq[-1]), origin='lower', aspect='auto', interpolation='nearest')
    p.colorbar(shrink=.5)
    p.xlabel('L')
    p.ylabel(r'$\nu$')

    p.subplot(224)
    plt_d = n.log10(n.abs(_d_orth)**2)
    mx = plt_d.max()
    p.imshow(plt_d, extent=(_L[0],_L[-1],_fq[0],_fq[-1]), vmax=mx, vmin=mx-drng, origin='lower', aspect='auto', interpolation='nearest')
    p.colorbar(shrink=.5)
    p.xlabel('u')
    p.ylabel(r'$\eta$')

p.show()




