#!/usr/bin/env python
#
#  print_fits_hdr.py
#  
#
#  Created by Danny Jacobs on 9/8/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse, pyfits as pf
from cStringIO import StringIO
o = optparse.OptionParser()

opts, args = o.parse_args(sys.argv[1:])

for file in args:
    hdulist = pf.open(file)
#    outfile = file[:-len('.fits')]+'.txt'
    for hdu in hdulist:
        outfile = StringIO()
        hdu.header.toTxtFile(outfile)
        outfile.reset()
        print ''.join(outfile.readlines())