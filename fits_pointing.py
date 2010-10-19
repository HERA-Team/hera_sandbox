#!/usr/bin/env python
#
#  fits_pointing.py
#  
#
#  Created by Danny Jacobs on 9/24/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pyfits as pf,ephem

o = optparse.OptionParser()
o.set_usage("fits_pointing.py *.fits")
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('-v',action='store_true',
    help='print more')
opts, args = o.parse_args(sys.argv[1:])

for file in args:
    hdulist = pf.open(file)
    RA = hdulist[0].header['CRVAL1']*n.pi/180
    DEC = hdulist[0].header['CRVAL2']*n.pi/180
    print ephem.Equatorial(RA,DEC).ra,ephem.Equatorial(RA,DEC).dec,
    if opts.v:print file
    else: print 