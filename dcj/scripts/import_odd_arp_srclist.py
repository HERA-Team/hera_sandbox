#!/usr/bin/env python
#
#  import_odd_arp_srclist.py
#  
#
#  Created by Danny Jacobs on 4/15/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])
print "#Name \t Ra \t Dec \t S_nu \t index \t mfreq \t Q"
for file in args:
    lines = open(file).readlines()
    for line in lines:
        if len(line.split(':'))>2: 
            ra = line.split('FLX')[0].split('_')[0].strip()
            dec = line.split('FLX')[0].split('_')[1].strip()[:-1]
        else:
            ra,dec = ('0','0')
        name = line.split('FLX')[0].strip()[:-1]
        S_nu = line.split('FLX=')[1].split('IND')[0].strip()
        index = line.split('IND=')[1].split('Q')[0]
        Q = line.split('Q=')[1].strip()
        print '\t'.join([name,ra,dec,S_nu,index,'0.15',Q])