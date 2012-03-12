#!/usr/bin/env python
#
#  trace.py
#  
#
#  Created by Danny Jacobs on 4/5/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--snap', dest='snap', action='store_true',
    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])