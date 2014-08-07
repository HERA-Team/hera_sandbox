#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import ephem, sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-m', '--map', dest='map',
    help='Haslam map.')
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='Channel')
opts,args = o.parse_args(sys.argv[1:])
fq = opts.freq

aa = a.cal.get_aa(opts.cal, .001, fq, 1)
im = a.img.Img(size=100, res=.5)
h = a.map.Map(fromfits=opts.map)

tx,ty,tz = im.get_top()
invalid = tx.mask.copy()
tx = tx.filled(0).flatten()
ty = ty.filled(0).flatten()
tz = tz.filled(0).flatten()

resp = aa[0].bm_response((tx,ty,tz), pol='y')**2
resp.shape = invalid.shape
resp = n.where(invalid, 0, resp/resp[0,0])

lsts = n.arange(0,2*n.pi,.01)
dat = []
for ha in lsts:
    print ha
    ex, ey, ez = im.get_eq(ra=ha, dec=aa.lat)
    ex = ex.filled(1).flatten() 
    ey = ey.filled(0).flatten()
    ez = ez.filled(0).flatten()
    tsky = h[ex,ey,ez] * (fq/.408)**-2.52
    tsky.shape = resp.shape
    tsky = n.where(invalid, 0, tsky)
    dat.append(n.sum(tsky * resp) / n.sum(resp))

dat = n.array(dat)
lsts = n.array(lsts)
lsts_fit = n.concatenate([lsts, lsts, lsts])
n.savez('tsys.npz',tsys=dat)
print dat
dat_fit = n.concatenate([dat, dat, dat])
sync_poly = n.polyfit(lsts_fit, dat_fit, deg=16)
print sync_poly
for ply in sync_poly: print float(ply)
sync_auto = n.polyval(sync_poly, lsts)

#p.subplot(211)
p.plot(lsts * 12 / n.pi, sync_auto)
p.plot(lsts * 12 / n.pi, dat)
p.show()

#resp = a.img.recenter(resp, (100,100))
#tsky = a.img.recenter(tsky, (100,100))

#print n.sum(tsky * resp) / n.sum(resp)

'''
p.subplot(131)
p.imshow(resp)
p.colorbar(shrink=.5)

p.subplot(132)
p.imshow(tsky)
p.colorbar(shrink=.5)

p.subplot(133)
p.imshow(tsky * resp)
p.colorbar(shrink=.5)
'''

p.show()
