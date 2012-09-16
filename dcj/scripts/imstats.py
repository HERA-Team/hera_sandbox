#!/usr/bin/env python
#
#  imstats.py
#  
#
#  Created by Danny Jacobs on 18 July 2012
#  PAPER Project
#


import aipy as a, numpy as n, math as m, ephem, sys,optparse, pyfits as pf
from kapteyn import wcs
from pylab import *
#import atpy
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('-r',dest='radius',default=10.,type='float',
    help="Analysis region around brightest source [default 10 deg]")
opts, args = o.parse_args(sys.argv[1:])
#import the image
#find the brightest source (print the location)
# find the dynamic range around this source

def query_disc(proj,pos,r):
    #input proj=wcs.Projection object for fits file
    #input pos = (RA,DEC) tuple in degrees
    #input r = search radius in degrees

    #the center pixel
    cpix = proj.topixel(pos)
    #the number of pixels to look at
    # = number spanned by the search area x pi
    pshape = (r/n.abs(proj.cdelt)*n.pi).astype(int)
    #form up the RA,DEC pixel positions in this area
    D,R = n.indices(pshape)
    S = n.sqrt((D - pshape[0]/2.)**2*n.abs(proj.cdelt[0])\
            + (R - pshape[1]/2.)**2*n.abs(proj.cdelt[1]))

    D += cpix[1] - pshape[1]/2.
    R += cpix[0] - pshape[0]/2.
    D.astype(int)
    R.astype(int)
    pix = n.array(zip(R[S<r],D[S<r]))
    #only return pixels that are within the image
    return pix[n.logical_and(pix[:,0]<proj.naxis[0],pix[:,1]<proj.naxis[1])]

aa = a.scripting.get_null_aa()
print "#filename pointingRA pointingDEC peakRA peakDEC peakFlux centralRMS DR"
for file in args:
    print "%s"%file,
    hdulist = pf.open(file)
    header = hdulist[0].header
    proj = wcs.Projection(header)
    image = n.ma.masked_invalid(hdulist[0].data.T)    
    try:
        map = proj.sub((proj.lonaxnum,proj.lataxnum))
    except:
        print "Aborting program. Could not find (valid) spatial map."
        raise
    maxpx = n.argwhere(image==image.max()).squeeze()
    F = image.max()
    RA,DEC = proj.toworld(maxpx)
    srcpos = n.array([RA,DEC])
    print proj.crval[0],proj.crval[1],RA,DEC,
    #get the pixels inside the surrounding annulus        
    center = n.array(image.shape).astype(int)/2
    doughpix = query_disc(proj,srcpos,opts.radius)
    #holepix = query_disc(proj,srcpos,1.)
    #donutpix = n.array([px for px in doughpix if px not in holepix])
    RMS = n.std(image[doughpix[:,1],doughpix[:,0]])   
    print F,RMS,F/RMS
