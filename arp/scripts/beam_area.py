#! /usr/bin/env python
"""
Uses *.hmap files to generate BeamAlm coefficients using the provided order
for a polynomial in frequency at each pointing, and the provided lmax and mmax
for fitting spherical harmonic coefficients to the spatial variation of each
of the polynomial coefficients.
"""

import aipy as a, numpy as n, sys, optparse, re

o = optparse.OptionParser()
o.set_usage('beam_area.py [options] *.hmap')
opts,args = o.parse_args(sys.argv[1:])
assert(len(args) > 0)

re_freq = re.compile(r'_(\d+)\.hmap$')
freqs = n.array([float(re_freq.search(f).groups()[0]) / 1e3 for f in args])
for i, filename in enumerate(args):
    print 'Reading', filename
    m = a.healpix.HealpixMap(fromfits=filename)
    m.map = (m.map / m[0,0,1]).clip(0,1)
    x,y,z = m.px2crd(n.arange(m.map.size), ncrd=3)
    v = n.where(z > 0, 1, 0)
    data = m.map.compress(v)
    print n.sum(data**2) * 4 * n.pi / m.npix()
    #x = x.compress(v)
    #y = y.compress(v)
    #z = z.compress(v)

