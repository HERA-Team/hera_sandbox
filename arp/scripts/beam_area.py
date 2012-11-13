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
o.add_option('-v', '--voltage_beam', action='store_true', help='Square data in file to produce power beams.')
opts,args = o.parse_args(sys.argv[1:])
assert(len(args) > 0)

re_freq = re.compile(r'_(\d+)\.hmap$')
#freqs = n.array([float(re_freq.search(f).groups()[0]) / 1e3 for f in args])
for i, filename in enumerate(args):
    print 'Reading', filename
    if True:
        m = a.healpix.HealpixMap(fromfits=filename)
    else:
        m2 = a.healpix.HealpixMap(fromfits=filename)
        m = a.healpix.HealpixMap(nside=32)
        m.from_hpm(m2)
        print m.nside()
    m.map = (m.map / m[0,0,1]).clip(0,1)
    x,y,z = m.px2crd(n.arange(m.map.size), ncrd=3)
    v = n.where(z > 0, 1, 0)
    beam = m.map.compress(v)
    if opts.voltage_beam: beam = beam**2
    nsum,nwgt = 0,0
    for i in range(100):
        noise = n.random.normal(size=beam.size)
        nsum += n.abs(n.sum(noise * beam))**2
        nwgt += 1
    #print n.sqrt(nsum / nwgt)
    print nsum / nwgt * 4 * n.pi / m.npix()
    #print beam.max(), beam.min(), float(v.sum()) / m.npix(), m.npix() / (4*n.pi)
    noise_beam = n.sum(beam) * 4 * n.pi / m.npix()
    #eor_beam = n.sum(beam**2)/n.sum(beam)**2 * m.npix() / (4*n.pi)
    eor_beam = n.sum(beam**2) * 4 * n.pi / m.npix()
    print u'   \u03A9 Nos:', noise_beam
    print u'   \u03A9 EoR:', eor_beam
    print u'   \u03A9 Eff:', noise_beam**2 / eor_beam
    
    #x = x.compress(v)
    #y = y.compress(v)
    #z = z.compress(v)

