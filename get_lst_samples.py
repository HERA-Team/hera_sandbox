#!/usr/bin/env python
#
#  get_lst_samples.py
#  
#
#  Created by Danny Jacobs on 6/2/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
import aipy as a, numpy as n
import sys

LST = 17*n.pi/12 #~2AM during May
DLST = 100 * a.const.arcmin
times = []
lsts = []
dat  = []
for f in sys.argv[1:]:
    print f,
    uv = a.miriad.UV(f)
    uv.select('antennae', 1, 10)
    #print "LST range",str(LST),"+/-",str(DLST)
    uv.select('ra',LST-DLST,LST+DLST)
    c = 0
    for (uvw,t,(i,j)),d in uv.all():
        times.append(t)
        lsts.append(uv['lst'])
        dat.append(d)
        c +=1
    print str(c)


times = n.array(times)
lsts = n.array(lsts)
dat = n.array(dat)
print times.shape,lsts.shape,dat.shape
dat = dat[:,500:550]
dat = n.average(dat,axis=1)
filename = "LST_%2.2f_%f_%f.dat" % (LST, n.min(times),n.max(times))

print times.shape,lsts.shape,dat.shape
#D = n.vstack((n.array(t),n.array(lsts),n.array(dat))).transpose()
D = n.vstack((times,lsts,dat)).transpose()
n.savetxt(filename,D)
