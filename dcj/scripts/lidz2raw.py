#!/usr/bin/env python
#
#  lidz2raw.py
#  
#
#  Created by Danny Jacobs on 8/5/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
from load_lidz import load_lidz
o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

for filename in args:
    F = load_lidz(filename)
    outname = filename[:-3]+'.raw'
    F.tofile(open(outname,'w'))
    del(F)