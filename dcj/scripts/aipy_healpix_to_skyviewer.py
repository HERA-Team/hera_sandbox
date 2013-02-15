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
#a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])


for file in args:
    outfile = '.'.join(file.split('.')[:-1])+'c.fits'
    print file + " > " + outfile
    hdulist = pf.open(file)
    hdu = hdulist[1]
    data = hdu.data.field('signal')[:]
    w = hdu.data.field('weights')[:]
    w = n.where(w>0,w,1)
    hdu.data.field('signal')[:] = data/w
    hdu.data.field('weights')[:] = n.ones_like(data)
    hdu.header['TTYPE1']='TEMPERATURE'
    hdu.header['TTYPE2']='N_OBS'
    hdu.header.update('COORDSYS','C')
    hdulist[1] = hdu

    hdulist.writeto(outfile,clobber=True)
