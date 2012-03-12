#!/usr/bin/env python
#
#  rfi_sum.py
#  
#
#  Created by Danny Jacobs on 12/16/09.
#  PAPER Project
#

import aipy as a, numpy as n, math as m
import sys, optparse,pickle as pkl

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])
for file in args:
    print file,
    jd = '.'.join(file.split('.')[0:2])
    mask = pkl.load(open(file))
    mask = n.array([mask[m] for m in mask])
    avg = n.average(mask,axis=0)
    outfile = jd+'.'+'avg'
    print " > ", outfile
    n.savetxt(outfile,avg,fmt='%2.3f')