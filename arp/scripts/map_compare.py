#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys

NSIDE = 512
m1 = a.map.Map(fromfits=sys.argv[-2])
m1.map.map *= 1.5
m1.to_fits('rescale.fits', clobber=True)
m2 = a.map.Map(fromfits=sys.argv[-1])

m2.map.map += m1.map.map
m2.wgt.map += m1.wgt.map
m2.to_fits('summap.fits', clobber=True)
sys.exit(0)

m1_m = a.healpix.HealpixMap(nside=NSIDE); m1_w = a.healpix.HealpixMap(nside=NSIDE)
m2_m = a.healpix.HealpixMap(nside=NSIDE); m2_w = a.healpix.HealpixMap(nside=NSIDE)
m1_m.from_hpm(m1.map); m1_w.from_hpm(m1.wgt)
m2_m.from_hpm(m2.map); m2_w.from_hpm(m2.wgt)
m1.map = m1_m; m1.wgt = m1_w
m2.map = m2_m; m2.wgt = m2_w

if True:
    valid = n.where(m1_w.map > m1_w.map.max() / 100, 1, 0)
    valid *= n.where(m2_w.map > m2_w.map.max() / 100, 1, 0)
    m1 = n.abs(m1_m.map / m1_w.map).compress(valid)
    m2 = n.abs(m2_m.map / m2_w.map).compress(valid)
    print n.sum(m1**2) / n.sum(m2**2)
    p.loglog(m1, m2, ',', alpha=.5)
    p.loglog([1e-6,1e6],[1e-6,1e6],'k-')
    p.show()

m1.to_fits('m1.fits', clobber=True)
m2.to_fits('m2.fits', clobber=True)
m1.map.map /= n.where(m2_m.map > 0, m2_m.map, 1)
m1.wgt.map /= n.where(m2_w.map > 0, m2_w.map, 1)
m1.to_fits('ratio.fits', clobber=True)
