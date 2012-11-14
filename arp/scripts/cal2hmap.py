#! /usr/bin/env python

import aipy as a, numpy as n
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])

FQ = .150

aa = a.cal.get_aa(opts.cal, n.array([FQ]))

h = a.healpix.HealpixMap(nside=64)
px = n.arange(h.npix())
x,y,z = h.px2crd(px)

resp = aa[0].bm_response((x,y,z))[0]
resp[0] = 1
import pylab; pylab.plot(resp); pylab.show()
print x[:10], y[:10], z[:10], resp[:10]
resp = n.where(z > 0, resp, 0)
h.set_map(resp) # writing a voltage beam
filename = opts.cal + '_%d.hmap' % (FQ*1e3)
print 'Writing', filename
h.to_fits(filename)
