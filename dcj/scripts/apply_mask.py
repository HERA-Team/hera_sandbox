#!/usr/bin/env python
#
#  apply_mask.py
#  
#
#  Created by Danny Jacobs on 2/3/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,os

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

def mfunc(uv,preamble,data,flags):
    return preamble, n.where(flags,0,data), flags

for file in args:
    uvofile = file+'M'
    print file +' > '+uvofile
    if os.path.exists(uvofile):
        print file, 'exists, skipping'
        continue
    uvi = a.miriad.UV(file)
    uvo = a.miriad.UV(uvofile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,raw=True,append2hist="APPLY_MASK: all masked values are set =0")
    