#!/usr/bin/env python
#
#  fix_healpy_fits.py
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


for file in args:
    outfile = '.'.join(file.split('.')[:-1])+'c.fits'
    print file + " > " + outfile
    hdulist = pf.open(file)
    hdu = hdulist[0]
    hdu.header.update('RESOLUTN', '9','Resolution index')
    hdu.header.update('TFORM1  ', '1E      ','4 bytes')
    hdulist[0] = hdu
    hdulist.writeto(outfile,clobber=True)