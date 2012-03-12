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
import pyfits as pf
o = optparse.OptionParser()
opts, args = o.parse_args(sys.argv[1:])

for file in args:
    hdulist = pf.open(file)
    print file
    print "\n".join(hdulist[0].header.get_history())
    print '-'*40