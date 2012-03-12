#!/usr/bin/env python
#
#  frankenfits.py
#  
#
#  Created by Danny Jacobs on 9/8/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pyfits as pf

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--snap', dest='snap', action='store_true',
    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

body = args[0]
brain = args[1]

bodyhdulist = pf.open(body)
braindulist = pf.open(brain)
bodyhdulist[0].header = braindulist[0].header
bodyhdulist[0].header['TFIELDS '] = 1
bodyhdulist[0].header['NAXIS '] = 1
bodyhdulist.writeto('arrrrg.fits',clobber=True)