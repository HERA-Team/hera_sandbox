#!/usr/bin/env python
#
#  fits_hrd2txt.py
#  
#
#  Created by Danny Jacobs on 9/8/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse, pyfits as pf
import cStringIO as StringIO
o = optparse.OptionParser()

opts, args = o.parse_args(sys.argv[1:])

for file in args:
    hdulist = pf.open(file)
#    outfile = file[:-len('.fits')]+'.txt'
    outfile = StringIO('')
    hdulist[0].header.toTxtFile(outfile)
    print ''.join(outfile.lines())