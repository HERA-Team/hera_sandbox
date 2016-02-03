#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C

def pk(k, k3pk=20.):
    return k3pk * 2*n.pi**2/k**3

def beam(sin_ang, sig=.3):
    return n.where(sin_ang**2 > 1, 0, n.exp(-sin_ang**2/(2*n.sin(sig)**2)))

z = 8.
kpl_mx = 0.5
bl_mx = 14e4
dk_du = C.pspec.dk_du(z)
dk_deta = C.pspec.dk_deta(z)
fq = C.pspec.z2f(z)
lam = a.const.c / fq / 1e9
print lam

bl_rng = n.arange(14e2, bl_mx, 14e2) # cm
kpr = (bl_rng/lam) * dk_du
kpr.shape = (1,kpr.size)
tau_rng = n.arange(0, kpl_mx/dk_deta, 10)
kpl = tau_rng * dk_deta
kpl.shape = (kpl.size,1)

ks = n.sqrt(kpl**2 + kpr**2)

pk_data = pk(ks)
#p.imshow(n.log10(pk_data.clip(1e4,1e8)),
#    extent=(bl_rng.min(),bl_rng.max(),0,kpl.max()),
#    origin='lower', interpolation='nearest', aspect='auto')
p.contour(bl_rng.flatten()/1e2, kpl.flatten(), n.log10(pk_data), vmax=15, vmin=3, linewidths=5)
#p.fill_between(bl_rng/1e2, bl_rng/a.const.len_ns*dk_deta, color='k', alpha=0.5)
p.plot(bl_rng/1e2, (bl_rng/a.const.len_ns+50)*dk_deta, 'k--')
p.plot(bl_rng/1e2, (bl_rng*n.sin(.14)/2/a.const.len_ns)*dk_deta, 'k:')
bl_rng.shape = (1,bl_rng.size)
sin_ang = kpl / dk_deta / (bl_rng / a.const.len_ns)
p.colorbar()
p.imshow(15.5+2*n.log10(beam(sin_ang)), vmax=15, vmin=3,
    extent=(bl_rng.min()/1e2,bl_rng.max()/1e2,0,kpl.max()),
    origin='lower', interpolation='nearest', aspect='auto')
p.colorbar()
p.xlim(14,300)
p.xlabel('Baseline length [m]')
p.ylim(0,0.5)
p.ylabel(r'$k_\parallel$', fontsize=16)
p.show()
