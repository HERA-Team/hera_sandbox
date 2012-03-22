#! /usr/bin/env python
import aipy as a, numpy as n
import sys

args = ['combmap_pgb011_psa010.fits', 'psa336_455_pgb399_smap_nice_v1c.fits']
maps = []
for filename in args:
    print 'Reading', filename
    m = a.map.Map(fromfits=filename)
    #_m = a.map.Map(nside=32)
    _m = a.map.Map(nside=16)
    _m.map.from_hpm(m.map)
    _m.wgt.from_hpm(m.wgt)
    maps.append(_m)

m1,m2 = maps
px = n.arange(m1.map.map.size)
x,y,z = m1.px2crd(px, ncrd=3)
rings = {}
for i in px:
    if z[i] > 0 and y[i] > 0: continue
    rings[z[i]] = rings.get(z[i], []) + [i]
zs, ratio = [], []
for z in rings:
    ring = n.array(rings[z])
    zs.append(z)
    m1val = n.where(m1[ring] > n.max(m1[ring])/10, 1, 0)
    m2val = n.where(m2[ring] > n.max(m2[ring])/10, 1, 0)
    
    #ratio.append(n.median(m1[ring] / m2[ring]).clip(.03,2))
    ratio.append(n.median(m1[ring].compress(m1val*m2val)) / n.median(m2[ring].compress(m1val*m2val)))

zs,ratio = n.array(zs), n.array(ratio)
valid = n.where(ratio > 0, 1, 0)
zs = zs.compress(valid)
ratio = ratio.compress(valid).clip(.1,10)

poly = n.polyfit(zs, n.log10(ratio), deg=12)
#poly = n.polyfit(n.abs(zs), n.log10(ratio), deg=6)
import pylab as p
p.plot(zs, ratio, '.')
p.plot(zs, 10**n.polyval(poly, zs), '.')
#p.plot(zs, 10**n.polyval(poly, n.abs(zs)), '.')
p.show()

m1 = a.map.Map(fromfits=args[0])
m1wgt = m1.wgt.map.max()
m1.map.map /= m1wgt; m1.wgt.map /= m1wgt
m2 = a.map.Map(fromfits=args[1])
m2wgt = m2.wgt.map.max()
m2.map.map /= m2wgt; m2.wgt.map /= m2wgt
px = n.arange(m1.map.map.size)
x,y,z = m.px2crd(px, ncrd=3)
gain = 10**n.polyval(poly, z)
#gain = 10**n.polyval(poly, n.abs(z))
#m1.map.map /= n.sqrt(gain)
#m2.map.map *= n.sqrt(gain)
m2.map.map *= gain
m2.map.map *= n.abs(1-y) * n.abs(z-.5)
m2.wgt.map *= n.abs(1-y) * n.abs(z-.5)
m1.map.map = m1.map.map + m2.map.map * gain
m1.wgt.map = m1.wgt.map + m2.wgt.map * gain
#prox = ((y-(-.1))/.1)**2 + ((z-0)/.2)**2
#m1.map.map = n.where(n.logical_and(y>-.1, z>0), m1.map.map,
#    n.where(z < -.2, m2.map.map, 
#    n.where(n.logical_and(y>-.2, z > -.2), (prox * m2.map.map + m1.map.map/gain**2),
#        m1.map.map/gain**2+m2.map.map)))
#m1.wgt.map = n.where(n.logical_and(y>-.1, z>0), m1.wgt.map, 
#    n.where(z < -.2, m2.wgt.map, 
#    n.where(n.logical_and(y>-.2, z > -.2), (prox * m2.wgt.map + m1.wgt.map/gain**2),
#        m1.wgt.map/gain**2+m2.wgt.map)))
m1.to_fits('rescale_'+filename, clobber=True)
    
