#!/usr/bin/env python
#
#  d_ogram.py
#  
#
#  Created by Danny Jacobs on 10/7/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

D = []
for file in args:
    uv = a.miriad.UV(file)
    aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
    for (uvw,t,(i,j)),dat in uv.all():
#        bl = a.miriad.ij2bl(i,j)
        uvw = aa.gen_uvw(i,j).squeeze()
        bl = n.average(uvw,axis=1)
        for d in dat:
            D.append((bl,d))
    