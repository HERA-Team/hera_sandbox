#!/usr/bin/env python
#
#  hpm_hpf.py
#  
#
#  Created by Danny Jacobs on 1/28/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse, healpy as hp
from scipy.special import erf,erfc
from pylab import *
o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('--lmin',dest='lmin',type='float',
    help='50% point for lower l limit.')
o.add_option('--lmax',dest='lmax',type='float',
   help='50% point for upper l limit.')
opts, args = o.parse_args(sys.argv[1:])

#generate the filter alms
#high-pass filter defined by:
#   error function with width of 50 l (hardcoded)
if not opts.lmin is None:
    dl_low = 40
    l0_low = opts.lmin
    filt_low  = lambda l: (erf((l-l0_low)*4/dl_low)+1)/2
    if opts.lmax is None: cl_filt = filt_low
if not opts.lmax is None:
    dl_hi = 50
    l0_hi = opts.lmax
    filt_hi = lambda l: erfc((l-l0_hi)*4/dl_hi)/2
    if opts.lmin is None: cl_filt = filt_hi
if not opts.lmin is None and not opts.lmax is None:
    cl_filt = lambda l: filt_hi(l)*filt_low(l)
l = n.linspace(opts.lmin-dl_low,opts.lmax+dl_low)
plot(l,cl_filt(l))
show()
for file in args:
    print file,'  >  ',
    s = hp.read_map(file)
    npix = len(s)
    nside = hp.npix2nside(npix)
    s_alm = hp.map2alm(s)
    l = n.linspace(0,3*nside-1)
    s_alm_filt = hp.almxfl(s_alm,cl_filt(l))
    s_filt = hp.alm2map(s_alm_filt,nside)
    outfile = file.split('/')[-1].split('.')
    if len(outfile)>2:
         outfile = '.'.join(outfile[:-1])+'f'+'.fits'
    else: outfile = '.'.join(outfile[:-1])+'.f.fits'
    print outfile
    hp.write_map(outfile,s_filt)