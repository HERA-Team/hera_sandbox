#! /usr/bin/env python

import numpy as n
import capo
import aipy
import glob
import sys, os
import optparse

""" 
Sorts UV files by LST order and returns the ordered list
"""

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

aa = aipy.cal.get_aa('psa6622_v003',n.array([0.15]))
files = args
lsts = []

for file in files:
    jd = float('.'.join(file.split('/')[-1].split('.')[1:3]))
    aa.set_jultime(jd)
    lst = aa.sidereal_time()
    lsts.append(lst)
order = n.argsort(lsts)
lsts_in_order = n.array(lsts)[order]
files_in_order =  n.array(files)[order]

for file in files_in_order:
    print file
