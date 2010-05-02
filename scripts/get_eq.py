#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, ephem as e, sys

data = [L for L in open(sys.argv[-1]).readlines() if not L.startswith('#')]
(lat, lon),data = data[0].split(), data[1:]
antpos1 = n.array([map(float, (L.split()+[0])[:3]) for L in data])
antpos1 = antpos1 - antpos1[0:1]
lat, lon = e.degrees(lat), e.degrees(lon)
print 'Using location:', lat, lon
m = a.coord.top2eq_m(0, lat)
print 'Topocentric coordinates:'
antpos1 /= a.const.len_ns / 1e2
print antpos1
print 'Equatorial coordinates:'
print n.dot(m, antpos1.transpose()).transpose()
x1,y1,z1 = antpos1[:,0], antpos1[:,1], antpos1[:,2]
p.plot(x1,y1, 'k.')

for i in range(len(x1)): p.text(x1[i], y1[i], str(i))
#p.xlim(-1000,0)
#p.ylim(-500,500)
p.show()

