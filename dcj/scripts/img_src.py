#!/usr/bin/env python
#
#  fits_history.py
#  
#
#  Created by Danny Jacobs on 9/7/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

