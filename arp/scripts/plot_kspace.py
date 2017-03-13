#! /usr/bin/env python

import numpy as n, pylab as p, capo as C, aipy as a

dk_deta = C.pspec.dk_deta(8.)
dk_du = C.pspec.dk_du(8.)

u = n.linspace(-1500,1500, 1024)
tau = u * 150. / a.const.len_ns

kpr = dk_du * u
kpl = dk_deta * tau


th = n.linspace(0,2*n.pi,1024)

p.plot(n.cos(th), n.sin(th), 'k', linewidth=2)
#p.plot(kpr, kpl, 'r', linewidth=2)
#p.plot(-kpr, kpl, 'r', linewidth=2)
#p.fill_betweenx(kpl, kpr, -kpr, alpha=.5, color='r', linewidth=1)
p.fill_between(kpr, kpl, -kpl, alpha=.5, color='r', linewidth=1)
p.xlim(-1,1)
p.xlabel(r'$k_\perp [h\ {\rm Mpc}^{-1}]$', fontsize=16)
p.ylabel(r'$k_\parallel[h\ {\rm Mpc}^{-1}]$', fontsize=16)
p.ylim(-1,1)
p.show()

