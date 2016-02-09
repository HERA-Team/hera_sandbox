import healpy as hp, numpy as n

#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p

#@p.ion()
#fqs = n.linspace(.1,.2,203)
fq = .15
bl1, bl2 = (0,26),(0,26)
N = 1   #number of universes to average over
aa = a.cal.get_aa('psa6240_v003', n.array([fq]))
h = a.healpix.HealpixMap(nside=64)

for i in xrange(N):
    print i
    sky = n.random.normal(size=h.map.size)
    h.map = sky # assume sky is in eq coord

fig = p.figure(1)
hp.mollview(h, min=-1, max=1, title='Unmasked map', fig=1, unit=r'$\Delta$T (mK)')
p.show()

