#! /usr/bin/env python
import aipy as a, numpy as n
import sys

for filename in sys.argv[1:]:
    h = a.map.Map(fromfits=filename)
    h.map.map /= h[0,0,1]
    px = n.arange(h.npix())
    x,y,z = h.px2crd(px, ncrd=3)
    data = h[x,y,z]**2 * n.where(z > 0, 1., 0)
    resp_px = data.sum()
    fwhm = n.where(data > .5, 1., 0).sum() * 4*n.pi/h.npix()
    #resp_px = n.where(z > 0, h[x,y,z]**4, 0).sum()
    #resp_px = n.where(n.logical_and(z > 0, h[x,y,z]**2 > .5), 1., 0).xum()
    #import pylab; pylab.plot(resp); pylab.show()
    px_size = 4*n.pi/h.npix()
    area = resp_px * px_size
    print filename, 'IntArea=%f FWHM=%f' % (area, fwhm)
