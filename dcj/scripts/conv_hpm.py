#!/usr/bin/env python
#
#  conv_hpm.py
#  
#
#  Created by Danny Jacobs on 1/30/10.
#  PAPER Project
#
"""
Convolve two healpix maps together.
conv_hpm.py map.fits psf.fits


"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
from healpy import *

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('--norm',dest='norm',action='store_true',
    help="normalize by the peak of the psf.")
opts, args = o.parse_args(sys.argv[1:])

map = read_map(args[0])
psf = read_map(args[1])
if opts.norm: psf = psf/n.max(psf)
psf_cl = alm2cl(map2alm(psf))
map_alm_fl = almxfl(map2alm(map),psf_cl)
map_f = alm2map(map_alm_fl,npix2nside(len(map)))
outfile = args[0].split('.')
if len(outfile)>2:
     outfile = '.'.join(outfile[:-1])+'f'+'.fits'
else: outfile = '.'.join(outfile[:-1])+'.f.fits'
print outfile
write_map(outfile,map_f)