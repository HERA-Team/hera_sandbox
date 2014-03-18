#! /usr/bin/env python
import numpy as n, pylab as p
import capo as C
import sys, re

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

re_z = re.compile(r'power_21cm_z(\d+\.\d+)\.dat')
npz = n.load('../fg_vs_umag_vs_fq.npz')

dat = {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    d = n.array([map(float, L.split()) for L in open(filename).readlines()])
    ks, pk = d[:,0], d[:,1]
    z = float(re_z.match(filename).groups()[0])
    k3pk = ks**3 / (2*n.pi**2) * pk# * mean_temp(z)**2
    dat[filename] = (ks, k3pk)

#zs = dat.keys(); zs.sort()
colors = 'kbgrcmy'
#fqs = n.arange(.150, .190, .01)
fqs = [.170]
for i,fq in enumerate(fqs):
    color = colors[i%len(colors)]
    for f in dat:
        p.loglog(dat[f][0], dat[f][1] * mean_temp(z)**2, color, label='%s,%f'%(f,z))

    # For 11x12 configuration
    #fq = .150
    B = .006
    NPOL = 2
    Trx = 150e3
    Tsky = 350e3 * (fq/.150)**-2.5
    Tsys = Trx + Tsky
    z = C.pspec.f2z(fq)
    print z, fq
    NDAYS = 120
    N_2HR_LSTBINS = 3.

    # Break up k modes sampled into parallel and perp components
    kpl_6MHz = C.pspec.dk_deta(C.pspec.f2z(fq)) * (1./B)
    kpl = 10**n.concatenate([n.arange(-20,-3,1), n.arange(-3, 1, .05)])
    kpr = C.pspec.dk_du(z) * 20 # assume 20-wavelength separation
    print 'k_pr', kpr
    ks = n.sqrt(kpl**2 + kpr**2)

    N_2HR_LSTBINS = 3.
    sense = 14.4551903993
    sense = sense * (ks / .1)**3 * (kpl / kpl_6MHz)**-0.5/ n.sqrt(N_2HR_LSTBINS) * (Tsys / 500e3)**2 * (120. / NDAYS) * (1. / NPOL)

    p.loglog(ks, sense, color+':')
#p.legend()
fq = .160
for ks, fqs, specs in zip(npz['ks'], npz['fqs'], npz['spec']):
    valid = n.where(n.abs(fqs - fq) < .002, 1, 0)
    valid *= n.where(specs == 1e6, 0, 1)
    print valid.shape, ks.shape, fqs.shape, specs.shape, valid.sum()
    valid = valid.flatten()
    ks = ks.flatten().compress(valid)
    specs = specs.flatten().compress(valid)
    print ks
    print specs
    p.loglog(n.abs(ks), n.abs(specs))
    p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', size=16)
    p.ylabel(r'$\Delta^2(k)\ [{\rm mK}^{2}]$', size=16)
    p.xlim(.06,4)
    p.ylim(1e-2, 1e5)
p.show()
