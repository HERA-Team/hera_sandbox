#!/usr/bin/env python

import numpy
import capo
import aipy
import pyfits
import pylab
import healpy

calfile = 'psa898_v003'
aa = aipy.cal.get_aa(calfile, numpy.array([0.15]))

gsm = '/Users/carinacheng/capo/ctc/images/gsm/gsm40/gsm1001.fits'

#view gsm map
img = healpy.read_map(gsm)
healpy.mollview(img, norm='log')
#pylab.show()

#get map coordinates
img = aipy.map.Map(fromfits=gsm, interp=True)
px = numpy.arange(img.npix()) #number of pixels
crd = numpy.array(img.px2crd(px,ncrd=3)) #coordinates for pixels
crds = numpy.asarray(crd) #equatorial
xcrd,ycrd,zcrd = crds[0],crds[1],crds[2]

#coordinate rotation
eq2top = aipy.coord.eq2top_m(aa.sidereal_time(),aa.lat)
t3 = numpy.dot(eq2top,crds)
tx,ty,tz = t3[0],t3[1],t3[2]

#get beam response
bmxx = aa[0].bm_response((tx,ty,tz),pol='x')**2
bmyy = aa[0].bm_response((tx,ty,tz),pol='y')**2
bm = numpy.where(tz > 0, 0.5 * (bmxx+bmyy),0)

#get image fluxes
imgfluxes = img[tx,ty,tz]

#plot gsm+beam
finalfluxes = numpy.reshape(imgfluxes*bm,(len(imgfluxes),))
wgts = numpy.ones_like(finalfluxes)
newmap = aipy.map.Map(nside=512)
newmap.add((xcrd,ycrd,zcrd),wgts,finalfluxes)
newmap.to_fits('beamimg.fits', clobber=True)

#fringe pattern
bl = aa.get_baseline(0,16)
blx,bly,blz = bl[0],bl[1],bl[2]
fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*.15)

#plot gsm+fringe
toplot = numpy.reshape(imgfluxes*fng*bm,(len(imgfluxes),))
newmap = aipy.map.Map(nside=512)
newmap.add((xcrd,ycrd,zcrd),wgts,toplot)
newmap.to_fits('fringeimg.fits', clobber=True)

