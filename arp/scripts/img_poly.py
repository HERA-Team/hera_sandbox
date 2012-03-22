#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, re

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, loc=True)
o.add_option('-d', '--deg', dest='deg', type='float', default=1,
    help='Degree of polynomial to use.')
o.add_option('--nogrid', dest='nogrid', action='store_true',
    help='Do not display RA/DEC grid.')

opts,args = o.parse_args(sys.argv[1:])

freqs = []
imgs = []
for filename in args:
    img,kwds = a.img.from_fits(filename)
    shape = img.shape[:2]
    print 'Reading %s' % filename
    freqs.append(kwds['freq'])
    print '   Freq:', freqs[-1]
    imgs.append(img.flatten().astype(n.float32))

freqs = n.array(freqs, dtype=n.float32)
imgs = n.array(imgs, dtype=n.float32)
print imgs.shape

print 'Fitting polynomials...'
print freqs.dtype, imgs.dtype
polys = n.polyfit(n.log10(freqs), n.log10(n.abs(imgs)), opts.deg)
print 'Done.'

print polys.shape
im0 = 10**n.polyval(polys, n.log10([.150])); im0.shape = shape
im1 = polys[0]; im1.shape = shape

import pylab as P

P.subplot(121)
if not opts.nogrid:
    from mpl_toolkits.basemap import Basemap
    xpx,ypx = shape
    dx1 = -(xpx/2 + .5) * kwds['d_ra'] * a.img.deg2rad
    dx2 = (xpx/2 - .5) * kwds['d_ra'] * a.img.deg2rad
    dy1 = -(ypx/2 + .5) * kwds['d_dec'] * a.img.deg2rad
    dy2 = (ypx/2 - .5) * kwds['d_dec'] * a.img.deg2rad
    map = Basemap(projection='ortho', lon_0=180, lat_0=kwds['dec'],
        rsphere=1, llcrnrx=dx1, llcrnry=dy1, urcrnrx=dx2,urcrnry=dy2)
    map.drawmeridians(n.arange(kwds['ra']-180,kwds['ra']+180,30))
    map.drawparallels(n.arange(-90,120,30))
    map.drawmapboundary()
    map.imshow(im0, interpolation='nearest')
else: P.imshow(im0, origin='lower', interpolation='nearest')
P.colorbar(shrink=.5)

P.subplot(122)
if not opts.nogrid:
    from mpl_toolkits.basemap import Basemap
    xpx,ypx = shape
    dx1 = -(xpx/2 + .5) * kwds['d_ra'] * a.img.deg2rad
    dx2 = (xpx/2 - .5) * kwds['d_ra'] * a.img.deg2rad
    dy1 = -(ypx/2 + .5) * kwds['d_dec'] * a.img.deg2rad
    dy2 = (ypx/2 - .5) * kwds['d_dec'] * a.img.deg2rad
    map = Basemap(projection='ortho', lon_0=180, lat_0=kwds['dec'],
        rsphere=1, llcrnrx=dx1, llcrnry=dy1, urcrnrx=dx2,urcrnry=dy2)
    map.drawmeridians(n.arange(kwds['ra']-180,kwds['ra']+180,30))
    map.drawparallels(n.arange(-90,120,30))
    map.drawmapboundary()
    map.imshow(im1, vmin=-3, vmax=1, interpolation='nearest')
else: P.imshow(im1, vmin=-3, vmax=1, origin='lower', interpolation='nearest')

P.colorbar(shrink=.5)

P.show()

#for i,filename in enumerate(args):
#    print '    ', filename
#    img,kwds = a.img.from_fits(filename)
#    im = simgs[i] ; im.shape = img.shape
#    a.img.to_fits(filename[:-8]+'sm'+filename[-8:], im, **kwds)
#    im = imgs[i] ; im.shape = img.shape
#    a.img.to_fits(filename[:-8]+'rs'+filename[-8:], im, **kwds)
#print 'Done.'
