#!/usr/bin/env python
#
#  ant_survey_stats.py
#  
#
#  Created by Danny Jacobs on 7/8/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#

import numpy as n,sys,aipy as a,optparse


o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=1)
o.add_option('--origin',default='a0 ',
    help='The entry (col 2 of csv trimble output file) to use as origin')
opts, args = o.parse_args(sys.argv[1:])

#self_cal=n.array([
#        [  8.79202085,      -1.6028028,     -12.23914682],
#        [294.264955417,    492.060738534,  -362.059217708],
#        [508.252844693,    355.993335315,  -623.206530202],
#        [487.575550872,   -316.824417864,  -611.383391464],
#
#        [ 38.6390644753,  -182.813356033,   -53.3660625273],
#        [  8.54237144392,   31.0509390628,  -11.3160595481],
#        [ 16.4643341327,   132.379151625,   -19.3325331207],
#        [120.66031028,     391.7692474,    -145.35850066],
#
#        [217.678227704,    472.106633788,  -264.597480588],
#        [440.047335423,    437.009263319,  -538.570719413],
#        [570.145149068,    199.378103243,  -704.823883029],
#        [561.376744633,   -157.265459422,  -702.063056458],
#
#        [402.782733597,   -397.899108667,  -512.388580721],
#        [269.923561549,   -431.303624138,  -350.450117519],
#        [191.417204318,   -403.457959543,  -250.462011768],
#        [ 87.5113892506,  -294.563815905,  -118.011949472],
#    ]).astype(n.float)
#
#self_cal -= n.array(self_cal[0])



#print "computing ",filename
aa=a.cal.get_aa(opts.cal, .1, .1, 1)

filename =args[0]
file = open(filename,'r')
data = [s.split(',') for s in file.readlines()]
ants = range(0,32)
#ants.remove(15)
#ants.remove(7)
#ants.remove(13)
ants = map(str,ants)
#ants.append('15-2')
#ants.append('7-2')
#ants.append('13-2')
#print ants
o = opts.origin
prefix = o[0]
p = n.zeros(3)
i = 0
names = []
for d in data:
        if d[0][0]!='#':
            if d[1].find(o)>0:
#                print d
                names.append(d[1][1:-1])
                p[0] += float(d[2]); p[1] += float(d[3])
                p[2] += float(d[4]) 
                i += 1
origin = p/i
print "origin found at:",names,origin
antpos = {}
pointcount = {}
names = []
#for ant in ants:
#    p = n.zeros(3)
#    i=0
for d in data:
    if d[0][0]!='#':
        name  = d[1][1:-1]
        if name.startswith(prefix):
            try:
                num = int(name[len(prefix):-1])
            except ValueError,e: 
                print "Input data error at %s. Point names must be <prefix><integer>"%\
                    d[1]
                print e
            if not d[1] in names: names.append(d[1])
            if not antpos.has_key(num): 
                antpos[num] = n.array(d[2:5]).astype(float)
                pointcount[num] = 1
            else: 
                antpos[num] += map(float,d[2:5])
                pointcount[num] += 1
print "found the following records:"
for name in names:
    print name
print antpos[1],pointcount[1]
nums = antpos.keys()
antpos = n.array([(antpos[ant]/pointcount[ant] - origin)*3.3356 for ant\
             in pointcount])
print "Raw ECEF coordinates"
print "{"
for cnt,p in enumerate(antpos):
#    print "%s,%3f,%3f,%3f,%d"%('a'+str(cnt),p[0],p[1],p[2],0)
     print "%d:[%f,%f,%f],"%(nums[cnt],p[0],p[1],p[2])
print "}"
long = aa.long
#long = -79.84997221765981*n.pi/180
m= a.coord.rot_m(long,n.array([0,0,1]))
antposr = n.dot(m,antpos.transpose()).transpose()
print "Raw ECEF coordinates rotated to",long," radians"
print "{"
for cnt,p in enumerate(antposr):
#    print "%s,%3f,%3f,%3f,%d"%('a'+str(cnt),p[0],p[1],p[2],0)
     print "%d:[%f,%f,%f],"%(nums[cnt],p[0],p[1],p[2])
print "}"
