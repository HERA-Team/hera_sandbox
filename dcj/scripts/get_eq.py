#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, ephem as e, sys

data = open(sys.argv[-1]).readlines()
(lat, lon),data = data[0].split(), data[1:]
antpos1 = n.array([map(float, (L.split()+[0])[1:4]) for L in data])
#print "absolute topocentric coordinates"
#print antpos1
antpos1 = antpos1 - antpos1[0:1]
lat, lon = e.degrees(lat), e.degrees(lon)
print 'Using location:', lat, lon
m = a.coord.top2eq_m(0, lat)
print 'Topocentric coordinates:'
antpos1[:,0] *= 674. / 656.
print antpos1
print 'Equatorial coordinates:'
print n.dot(m, antpos1.transpose()).transpose()
x1,y1,z1 = antpos1[:,0], antpos1[:,1], antpos1[:,2]
#p.plot(x1,y1, 'k.')

data = open(sys.argv[-1]).readlines()
(lat, lon),data = data[0].split(), data[1:]
antpos2 = n.array([map(float, (L.split()+[0])[0:3]) for L in data])
antpos2 = antpos2 - antpos2[0:1]
lat, lon = e.degrees(lat), e.degrees(lon)
m = a.coord.eq2top_m(0, lat)
antpos2 = n.dot(m, antpos2.transpose()).transpose()
x2,y2,z2 = antpos2[:,0], antpos2[:,1], antpos2[:,2]
#p.plot(x2,y2, 'r.')
#for i in range(len(x1)): p.text(x1[i], y1[i], str(i))
#p.xlim(-1000,0)
#p.ylim(-500,500)
#p.show()

