#!/usr/bin/env python
#
#  rfi_avg_compile.py
#  
#
#  Created by Danny Jacobs on 12/16/09.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pickle as pkl

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

def jd2hr(jd):
    return (jd % int(jd))*24+12

occupancy_jd = {}
hours = {}
for file in args:
    print file
    imp = n.loadtxt(file)
    jd = float('.'.join(file.split('.')[0:2]))
    occupancy_jd[jd] = imp
    hours[jd] = jd2hr(jd)
print "  > ","rfi_occupancy_jd.txt"
#times = n.sort(occupancy_jd.keys())
#occupancy_jd = n.array([occupancy_jd[t] for t in times])
#n.savetxt('rfi_occupancy_jd.txt',occupancy_jd)
file = open('rfi_occupancy_jd.txt','w')
pkl.dump(occupancy_jd,file)    