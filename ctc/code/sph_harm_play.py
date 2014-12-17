#!/usr/bin/env python

"""

NAME: 
      sph_harm_play.py 
PURPOSE:
      -Play around with spherical harmonics & coefficients
EXAMPLE CALL:
      ./sph_harm_play.py
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt
import optparse
from scipy import special
import os, sys
import healpy

LMAX = 500 #highest is 4*nside?
MMAX = 500

#playing with Healpy

#img = aipy.map.Map(nside=512)
#img.map.map = numpy.zeros_like(img.map.map)
#img.put((0,0,1),1.0,1.0) #use eq coordinates
#img.to_fits('/Users/carinacheng/capo/ctc/images/test2.fits', clobber=True)
#img = healpy.read_map('/Users/carinacheng/capo/ctc/images/test2.fits')

img = healpy.read_map('/Users/carinacheng/capo/ctc/images/gsm/gsm20/gsm1001.fits')

#healpy.mollview(img, norm='log')
alms = healpy.sphtfunc.map2alm(img,lmax=LMAX,mmax=MMAX)
print alms
print len(alms)
#print alms
alm_img = healpy.sphtfunc.alm2map(alms,nside=512)#,lmax=LMAX,mmax=MMAX)
#print alm_img
#healpy.mollview(alm_img, norm='log')

cls = healpy.sphtfunc.alm2cl(alms)
#plt.plot(cls)

#pylab.show()



"""
#plotting alm's of a Healpix Map

#img = aipy.map.Map(nside=512)
#img.map.map = numpy.zeros_like(img.map.map)
#img.put((0,0,1),1.0,100.0) #use eq coordinates

img = aipy.map.Map(fromfits = '/Users/carinacheng/capo/ctc/images/gsm/gsm20/gsm1001.fits')

alms = img.to_alm(lmax,mmax) #(lmax, mmax)
#print alm[1,0]
img.from_alm(alms)
img.to_fits('/Users/carinacheng/capo/ctc/images/test2.fits', clobber=True)
alm_img = healpy.read_map('/Users/carinacheng/capo/ctc/images/test2.fits')
healpy.mollview(alm_img)
"""

"""
#plotting spherical harmonics on Healpix Map

img = aipy.map.Map(nside=128)
px = numpy.arange(img.npix())
thetas,phis = img.px2crd(px,ncrd=2)

for i in range(len(thetas)):
    
    value = special.sph_harm(2,2,thetas[i],phis[i]) #(mmax, lmax, theta_crd, phi_crd)
    img.put((thetas[i],phis[i]),1.0,value.real)

img.to_fits('/Users/carinacheng/capo/ctc/images/test2.fits', clobber=True)
"""
