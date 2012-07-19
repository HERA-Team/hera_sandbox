#!/usr/bin/env python
#
#  fix_aipy_fits.py
#  
#
#  Created by Danny Jacobs on 18 July 2012
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pyfits as pf

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])


for file in args:
    outfile = '.'.join(file.split('.')[:-1])+'c.fits'
    print file + " > " + outfile
    hdulist = pf.open(file)
    hdu = hdulist[0]
    hdu.header['CTYPE1']='RA---SIN'
    hdulist[0] = hdu

    hdulist.writeto(outfile,clobber=True)
