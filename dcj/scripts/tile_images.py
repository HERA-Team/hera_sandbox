#! /usr/bin/env python
from astropy.io import fits
import numpy as n
import sys

for filename in sys.argv[1:]:
    tilefile = filename[:-4]+'_tile.fits'
    print filename+'>'+tilefile
    hdulist = fits.open(filename)
    data = hdulist[0].data
    bigdata = n.tile(data,(1,3,3))
    hdulist[0].data = bigdata
    hdulist[0].header['CRPIX1'] = int(bigdata.shape[2]/2.)
    hdulist[0].header['CRPIX2'] = int(bigdata.shape[1]/2.)
    hdulist.writeto(tilefile,clobber=True)
