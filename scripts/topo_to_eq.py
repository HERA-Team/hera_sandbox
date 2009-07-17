#!/usr/bin/env python
#
#  topo_to_eq.py
#  
#
#  Created by Danny Jacobs on 7/12/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#

import aipy as a, numpy as n, pylab as p, ephem as e, sys,optparse


o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--csv',dest='csv',action='store_true',
    help='Print csv instead of cal file output')
opts, args = o.parse_args(sys.argv[1:])

data = open(args[0]).readlines()
aa=a.cal.get_aa(opts.cal, .1, .1, 1)
data= [L.split(',')[1:4] for L in data]

antpos1 = n.array([map(float, L) for L in data])
antpos1 = (antpos1 - antpos1[0])*3.3356#ns/m
antpos1 = n.array([[pos[1],pos[0],pos[2]] for pos in antpos1]).squeeze()
lat, lon = e.degrees(aa.lat), e.degrees(aa.long)
if not opts.csv:
    print 'Using location:', lat, lon
    m = a.coord.top2eq_m(0, lat)
    print 'Topocentric coordinates:'
    print antpos1
    print 'Equatorial coordinates:'
    eq= n.dot(m, antpos1.transpose()).transpose()
    print "[",
    for q in eq:
        print "[%f,%f,%f],"%(q[0],q[1],q[2])
    print "]"
else:
    m = a.coord.top2eq_m(0, lat)
    eq= n.dot(m, antpos1.transpose()).transpose()
    for i,q in enumerate(eq):
        print "a%d,%f,%f,%f"%(i,q[0],q[1],q[2])
#x,y,z = antpos1[:,0], antpos1[:,1], antpos1[:,2]
#print zip(x,y,z)
#p.plot(x1,y1, 'k.')
#
#data = open(sys.argv[-1]).readlines()
#(lat, lon),data = data[0].split(), data[1:]
#antpos2 = n.array([map(float, (L.split()+[0])[0:3]) for L in data])
#antpos2 = antpos2 - antpos2[0:1]
#lat, lon = e.degrees(lat), e.degrees(lon)
#m = a.coord.eq2top_m(0, lat)
#antpos2 = n.dot(m, antpos2.transpose()).transpose()
#x2,y2,z2 = antpos2[:,0], antpos2[:,1], antpos2[:,2]
#p.plot(x2,y2, 'r.')
#for i in range(len(x1)): p.text(x1[i], y1[i], str(i))
#p.xlim(-1000,0)
#p.ylim(-500,500)
p.show()
