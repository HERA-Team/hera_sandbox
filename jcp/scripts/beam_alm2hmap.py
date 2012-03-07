#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))
im = a.img.Img(size=100, res=.5)

tx,ty,tz = im.get_top()
invalid = tx.mask.copy()
tx = tx.filled(0).flatten()
ty = ty.filled(0).flatten()
tz = tz.filled(0).flatten()

resp = aa[0].bm_response((tx,ty,tz), pol='y')**2
resp.shape = invalid.shape
resp = n.where(invalid, 0, resp/resp[0,0])
resp.shape = tx.shape

beam = a.map.Map(nside=64,interp=True)
beam.add((tx,ty,tz),1,resp)
beam.to_fits('alm_model150mhz_n64.fits',clobber=True)
