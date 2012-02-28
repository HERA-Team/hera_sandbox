#!/usr/bin/python

import aipy, readtle, numpy as n
import sys, math
import pylab as p
import optparse, ephem

o = optparse.OptionParser()
opts, args = o.parse_args(sys.argv[1:])

aa = readtle.get_aa(n.array([.1375]))
aa.set_jultime(2454952.88)
#cat = readtle.get_catalog(['orb_22'])
cat = readtle.get_catalog()

cat.compute(aa)
src = []
for s in cat:
    src.append(cat[s])

az = []
alt = []
ra = []
dec = []
sublong = []
sublat = []

for i in range(len(src)):
    az.append([])
    alt.append([])
    ra.append([])
    dec.append([])
    sublat.append([])
    sublong.append([])

lat = ephem.degrees('38:25:59.24')
lon = ephem.degrees('-79:51:02.1')

decimator = 0
for uvfile in args:
    print 'Reading', uvfile
    #uv = aipy.miriad.UV(sys.argv[-1])
    uv = aipy.miriad.UV(uvfile)
    uv.select('auto',0,0)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if decimator == 0:
            decimator += 1
            aa.set_jultime(t)
            cat.compute(aa)
            for c,r in enumerate(cat):
                #ra[c].append(cat[r].ra)
                #dec[c].append(cat[r].dec)
                ra[c].append(cat[r].a_ra)
                dec[c].append(cat[r].a_dec)
                sublat[c].append(cat[r].sublat)
                sublong[c].append(cat[r].sublong)
                alt[c].append(cat[r].alt)
                az[c].append(cat[r].az)
        decimator += 1
        decimator = decimator % 1000

def degrees(a): return math.degrees(a)
p.subplot(311)
for i in range(len(src)):
    p.plot(map(degrees,sublong[i]),map(degrees,sublat[i]),'.')
p.plot([math.degrees(lon)],[lat],'o')
p.subplot(312)
for i in range(len(src)):
    p.plot(map(degrees,ra[i]),dec[i],'.')
p.subplot(313)
for i in range(len(src)):
    p.plot(map(degrees,az[i]),map(degrees,alt[i]),'.')
p.show()

"""p.subplot(111, polar=True)
test_az = n.arange(360)
test_alt = n.ones_like(test_az) * 0.1
p.plot(test_az,test_alt)
test_alt = n.ones_like(test_az) * 0.5
p.plot(test_az,test_alt)
test_alt = n.ones_like(test_az) * 0.01
p.plot(test_az,test_alt)
p.show()"""

# def rt_2_x(r,t): return float(r) * math.cos(abs(float(t)))
# def rt_2_y(r,t): return float(r) * math.sin(abs(float(t)))
# def xy_2_r(x,y): return math.sqrt(x*x + y*y)
# def xy_2_t(x,y): return math.atan(y/x)
# def fix_alt(alt): return float(repr(alt))
# def fix_az(az): return float(repr(az))*360.0/(2*math.pi)
# 
# #p.subplot(211)
# for i in range(len(src)):
#     #alt2 = map(fix_alt,alt[i])
#     #az2 = map(fix_az,az[i])
#     p.subplot(111, polar=True)
#     p.plot(alt[i],az[i])
# #    p.plot(range(len(alt[i])),alt[i])
#      #x = map(rt_2_x, alt[i],az[i])
#      #y = map(rt_2_y, alt[i],az[i])
#      #p.plot(x,y)
#      #p.plot(az[i],alt2)
# #p.subplot(212)
# #for i in range(len(src)):
# #    p.plot(range(len(az[i])),az[i])
# #ax = p.subplot(111, polar=True)
# #for i in range(len(src)):
# #    #r = math.sqrt(math.pow(az,2) + math.pow(alt,2))
# #    #theta = math.atan(math.abs(alt), math.abs(az))
# #    p.plot(theta, r)
# p.show()

