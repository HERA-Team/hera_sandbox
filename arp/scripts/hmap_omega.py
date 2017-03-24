#! /usr/bin/env python
import aipy as a, numpy as n
import sys

for filename in sys.argv[1:]:
    h = a.map.Map(fromfits=filename)
    print h.nside()
    h.map.map /= h[0,0,1]
    aa = a.cal.get_aa('psa6240_v003', n.array([.150]))
    px = n.arange(h.npix())
    x,y,z = h.px2crd(px, ncrd=3)
    #data = h[x,y,z]**2 * n.where(z > 0, 1., 0)
    data = aa[0].bm_response((x,y,z),pol='x')[0]**2 * n.where(z > 0, 1., 0)
    resp_px = data.sum()
    resp2_px = n.sum(data**2)
    fwhm = n.where(data > .5, 1., 0).sum() * 4*n.pi/h.npix()
    #resp_px = n.where(z > 0, h[x,y,z]**4, 0).sum()
    #resp_px = n.where(n.logical_and(z > 0, h[x,y,z]**2 > .5), 1., 0).xum()
    #import pylab; pylab.plot(resp); pylab.show()
    px_size = 4*n.pi/h.npix()
    area = resp_px * px_size
    area2 = resp2_px * px_size
    print filename, 'IntArea=%f IntArea^2=%f FWHM=%f' % (area, area2, fwhm)
    print area**2 / area2
