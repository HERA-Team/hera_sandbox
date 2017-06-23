#! /usr/bin/env python
"""
Uses *.hmap files to generate BeamAlm coefficients using the provided order
for a polynomial in frequency at each pointing, and the provided lmax and mmax
for fitting spherical harmonic coefficients to the spatial variation of each
of the polynomial coefficients.
"""

import aipy as a, numpy as n, pylab as p, sys, optparse, re
from mpl_toolkits.basemap import Basemap
from matplotlib import rc
import healpy as hp
rc('text', usetex=True)

o = optparse.OptionParser()
o.set_usage('hmap2alms.py [options] *.hmap')
a.scripting.add_standard_options(o, pol=True)
o.set_description(__doc__)
o.add_option('--order', dest='order', type=int, default=7,
    help='Polynomial order to fit for each pixel vs. frequency.  Default 7')
o.add_option('--lmax', dest='lmax', type=int, default=8,
    help='LMAX parameter for fitting spherical harmonics. Default 8')
o.add_option('--mmax', dest='mmax', type=int, default=8,
    help='MMAX parameter for fitting spherical harmonics. Default 8')
o.add_option('--abssq', action='store_true',
    help='Takes the absolute value squared of the input beams. Appropriate for making power beams.')
o.add_option('--mirror', action='store_true',
    help='Enforces mirror symmetry of the beams.')
opts,args = o.parse_args(sys.argv[1:])
print sys.argv[1:]
assert(len(args) > 0)

nplt = len(args)
RES = .01
SZ = (int(n.pi/2/RES), int(2*n.pi/RES)+2)
alt,az = n.indices(SZ)
alt = n.pi/2 - alt.astype(n.float) * RES
az = az.astype(n.float) * RES
px,py,pz = a.coord.azalt2top((az.flatten(), alt.flatten()))
bmap = Basemap(projection='ortho', lat_0=90, lon_0=180)
cx,cy = bmap(az * a.img.rad2deg, alt * a.img.rad2deg)

def plot_hmap(hmap, cnt, max=0, min=-2):
    p.subplot(3, nplt, cnt)
    hmap.set_interpol(True)
    data = hmap[px,py,pz]
    data = n.log10(data**2)
    data.shape = az.shape
    bmap.drawmapboundary()
    bmap.drawmeridians(n.arange(0,360,30))
    bmap.drawparallels(n.arange(0,90,10))
    step = (max - min) / 35.
    levels = n.arange(min-step, max+step, step)
    bmap.contourf(cx, cy, data, levels)

re_freq = re.compile(r'_(\d+)\.hmap$')
freqs = n.array([float(re_freq.search(f).groups()[0]) / 1e3 for f in args])
hmaps = []
dat = []
for i, filename in enumerate(args):
    print 'Reading', filename
    m = a.healpix.HealpixMap(fromfits=filename)
    if opts.abssq: m.map = n.abs(m.map)**2
    m.map = (m.map / m[0,0,1]).clip(0,1)
    if len(dat) == 0:
        x,y,z = m.px2crd(n.arange(m.map.size), ncrd=3)
        v = n.where(z > 0, 1, 0)
        x = x.compress(v)
        y = y.compress(v)
        z = z.compress(v)
    # Swap pol from file 
    if opts.pol == 'y': m[y,x,z] = m[x,y,z]
    if opts.mirror:
        m[x,y,z] = m[x,y,z]/2.0 + m[-x,y,z]/2.0
        m[x,y,z] = m[x,y,z]/2.0 + m[x,-y,z]/2.0
    # Mirror bottom to prevent sharp discontinuities at z=0
    m[x,y,-z] = m[x,y,z]
    plot_hmap(m, i+1)
    p.title('%0.3f GHz' % freqs[i])
    
    hmaps.append(m)
    dat.append(m.map)
    nside = m.nside()
dat = n.array(dat)

print 'Fitting polynomials in frequency at each pointing...'
poly = []
for i in range(dat.shape[1]):
    poly.append(n.polyfit(freqs, dat[:,i], opts.order))
poly = n.array(poly)

print 'Creating Healpix maps of each polynomial coefficient'
pmaps = []
for i in range(opts.order+1):
    pmap = a.healpix.HealpixMap(nside, scheme='RING')
    pmap.map = poly[:,i]
    pmap.set_interpol(True)
    pmaps.append(pmap)

print 'Computing spherical harmonic coefficients...'
palms = [pmap.to_alm(opts.lmax, opts.mmax) for pmap in pmaps]
    
print 'Smoothing polynomial data...'
for i, pmap in enumerate(pmaps):
    print 'ALM Coeffs for Poly Term', opts.order - i
    data = palms[i].get_data()
    data = n.transpose(n.array([data.real, data.imag])).flatten()
    print list(data)
    print '-----------------------------------------------------------'
    pmap.from_alm(palms[i])

dat = n.array([pm.map for pm in pmaps])
resp = []
for i in range(dat.shape[1]):
    resp.append(n.polyval(dat[:,i], freqs))
resp = n.array(resp)
for i in range(len(args)):
    m.map = resp[:,i].clip(0,1)
    plot_hmap(m, len(args)+i+1)
    hmaps[i].map = n.abs(hmaps[i].map - m.map)
    plot_hmap(hmaps[i], 2*len(args)+i+1, max=-2, min=-4)

p.show()

