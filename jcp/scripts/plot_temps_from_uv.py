#! /usr/bin/env python

import numpy as n, aipy as a, pylab as p, optparse, sys

o = optparse.OptionParser()
opts, args = o.parse_args(sys.argv[1:])

ts = []
t_cable, t_balun, t_load = [],[],[]
for file in args:
    try: uv = a.miriad.UV(file)
    except(IOError): continue
    print file
    uv.select('antennae',0,1)
    for (crd,t,(i,j)), d in uv.all():
        ts.append(t)
        t_cable.append(uv['t_cable'])
        t_balun.append(uv['t_balun'])
        t_load.append(uv['t_load'])

p.subplot(131)
p.plot(ts,t_cable)
p.subplot(132)
p.plot(ts,t_balun)
p.subplot(133)
p.plot(ts,t_load)
p.show()

