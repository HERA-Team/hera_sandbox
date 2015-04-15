#! /usr/bin/env python

import numpy as n, pylab as p, capo as C
npz = n.load('data/coverage.npz')
lsts,cnt,var = npz['times'], npz['cnt'], npz['var']
lsts = n.where(lsts > 5, lsts - 2*n.pi, lsts) * 12 / n.pi
fqs = n.linspace(.1,.2,var.shape[1])
jy2T = C.pspec.jy2T(fqs)
jy2T.shape = (1,203)
extent = (1e3*fqs[0], 1e3*fqs[-1], lsts[-1], lsts[0])
Tlst = n.sqrt(var)*jy2T
#three_sig_correction = 0.974
three_sig_correction = 1/1.34
Tsys = Tlst * 1e-3 * n.sqrt(43 * 100. / 203 * 1e6) / three_sig_correction
Trms = Tlst / n.sqrt(cnt.clip(1,n.Inf))

p.figure(figsize=(7,4.2))
chunk = 83
for i in xrange(40,600,chunk):
    print i, n.average(lsts[i:i+chunk])
    lst_hr = int(n.around(n.average(lsts[i:i+chunk])))
    Tsys_avg = n.sqrt(n.average(Tsys[i:i+chunk]**2, axis=0))
    valid = n.where(Tsys_avg > 213, 1, 0)
    p.plot(1e3*fqs.compress(valid), Tsys_avg.compress(valid), label='LST = %02dh00'%lst_hr)
p.legend(loc='best', fontsize='medium')
p.xlim(125,175)
p.ylim(400,800)
p.grid()
p.xlabel(r'${\rm Frequency\ [MHz]}$', fontsize=16)
p.ylabel(r'${\rm T}_{\rm sys}$', fontsize=16)
p.subplots_adjust(0.12, 0.15, .95, .95)
p.show()

Tsys_jy = Tsys / 1e-3 / jy2T

#n.savez('tsys_model.npz', tsys=Tsys, freq=fqs, lsts=lsts)
#n.savez('tsys_model_jy.npz', tsys_jy=Tsys_jy, freq=fqs, lsts=lsts)
if True:
    p.subplot(121)
    C.arp.waterfall(Tsys, mx=700, drng=500, mode='lin', extent=extent)
    p.xlim(110,185)
    p.xlabel('Frequency [MHz]')
    p.xticks([125, 150, 175])
    p.colorbar()

    p.subplot(122)
    C.arp.waterfall(Trms, mx=30, mode='lin', extent=extent)
    p.xlim(110,185)
    p.xlabel('Frequency [MHz]')
    p.xticks([125, 150, 175])
    p.colorbar()
    p.show()
