#!/usr/bin/env python
#
#  sum_fits.py
#  
#
#  Created by Danny Jacobs on 6/24/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse, pyfits,ephem
"""
Sum input fits files.
sum_fits.py *.fits -o sum.fits
"""
o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
o.add_option('--avg',action='store_true',
    help='Average instead of summing.')
o.add_option('-o', dest='outfile', default='newstack.fits',
    help='Output file.')
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

for i,slice in enumerate(args):
    print "loading: "+slice
    data, kwds = a.img.from_fits(slice)
    if i==0:
        outdata = data
    else:
        outdata += data
if opts.avg: outdata /= i

files = ["%s  %6.2f  %6.2f"%(l.split('/')[-1],kwds['ra'],kwds['dec']) for l in args]
history = "stack_fits.py:"+ '\n  '+'\n  '.join(files)
a.img.from_fits_to_fits(args[0],opts.outfile,outdata,kwds,
    history=history)
print 'created: '+opts.outfile
