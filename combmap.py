#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-j','--juldates', dest='juldates', 
    help='A list of Julian dates corresponding to when each map was observed.  Used to downweight a region around the Sun in each map.  Entering 0 for JD omits Sun downweighting for that map')
o.add_option('-s','--sunavoid', dest='sunavoid', type='float', default=20.,
    help='Radius (in degrees) of gaussian down-weighting around location of Sun. Default 20')
opts,args = o.parse_args(sys.argv[1:])

sun = a.src.get_catalog(['Sun'])
sun = sun.values()[0]
juldates = map(float, opts.juldates.split(','))
sunavoid = opts.sunavoid * a.img.deg2rad
aas = [a.cal.get_aa(c, .001, .150, 1) for c in opts.cal.split(',')]

m = None
for aa, map, j in zip(aas, args, juldates):
    print 'Reading', map
    _m = a.map.Map(fromfits=map)
    scale = _m.wgt.map.sum()
    th,phi = _m.px2crd(n.arange(_m.npix()), ncrd=2)
    dec = n.pi/2 - th
    alt = n.pi/2 - n.abs(dec - aa.lat)
    az = n.zeros_like(alt)
    tx,ty,tz = a.coord.azalt2top((az,alt))
    resp = aa[0].beam.response((tx,ty,tz)).squeeze()
    print aa.lat * a.img.rad2deg, resp.shape
    resp = n.where(alt > 0, resp**2, 0)
    resp /= resp.sum()
    if j == 0: swgt = 1
    else:
        aa.set_jultime(j)
        sun.compute(aa)
        print sun.ra, sun.dec, j
        ex,ey,ez = _m.px2crd(n.arange(_m.npix()), ncrd=3)
        sx,sy,sz = sun.get_crds('eq', ncrd=3)
        r2 = (ex-sx)**2 + (ey-sy)**2 + (ez-sz)**2
        swgt = n.where(r2 < 1, 1-n.exp(-r2 / sunavoid**2), 1)
    N = len(aa)
    resp *= N * (N+1) / 2
    _m.map.map *= resp * swgt / scale
    _m.wgt.map *= resp * swgt / scale
    #_m.map.map = resp * swgt
    #_m.wgt.map = n.ones_like(resp)
    if m is None: m = _m
    else:
        m.map.map += _m.map.map
        m.wgt.map += _m.wgt.map
m.to_fits('combmap.fits', clobber=True)
