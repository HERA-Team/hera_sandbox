#!/usr/bin/env python
#
#  img_cat_flux.py
#  
#
#  Created by Danny Jacobs on 10 April 2012.
#  PAPER Project
#


import aipy as a, numpy as n, math as m, ephem, sys,optparse, pyfits as pf
from kapteyn import wcs
from pylab import *
import atpy

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('-r',dest='radius',default=1.,type='float',
    help="Analysis region around each input source, in degrees. [1]")
o.add_option('-o',dest='outfile',
    help="Output the catalog in a 'standard format' [None]")

opts, args = o.parse_args(sys.argv[1:])

outer_radius = 4
inner_radius = 1
n.set_printoptions(precision=2,suppress=True)

def pixsep(proj,pxA,pxB):
    A = n.array(proj.toworld(pxA)).astype(float)
    B = n.array(proj.toworld(pxB)).astype(float)
    return B-A

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
    S = n.sqrt((D - pshape[0]/2.)**2*proj.cdelt[0]\
            + (R - pshape[1]/2.)**2*proj.cdelt[1])

    D += cpix[1] - pshape[1]/2.
    R += cpix[0] - pshape[0]/2.
    D.astype(int)
    R.astype(int)

#    #compute the RA/DEC distances from center
#    dD,dR = n.zeros_like(D).astype(float),n.zeros_like(R).astype(float)
#    for ra in range(pshape[0]):
#        for dec in range(pshape[1]):
#            try:
#                dR[ra,dec],dD[ra,dec] = pixsep(proj,cpix,
#                    (R[ra,dec],D[ra,dec]))
#            except(wcs.WCSinvalid):
#                #point is not on the image
#                dR[ra,dec],dD[ra,dec] = 180,180 #put it as far away as possible
#    S = n.sqrt(dR**2 + dD**2)
    pix = n.array(zip(R[S<r],D[S<r]))
    #only return pixels that are within the image
    print len(n.logical_and(pix[:,0]<proj.naxis[0],pix[:,1]<proj.naxis[1]))
    return pix[n.logical_and(pix[:,0]<proj.naxis[0],pix[:,1]<proj.naxis[1])]
    #return S,R,D

clist,ccoff,ccats = a.scripting.parse_srcs(opts.src,opts.cat)
print clist,ccoff,ccats
if not opts.cal is None:
    print "cal"
    ccat = a.cal.get_catalog(opts.cal,srcs=clist,cutoff=ccoff,catalogs=ccats)
else:
    ccat = a.src.get_catalog(srcs=clist,cutoff=ccoff,catalogs=ccats)
aa = a.scripting.get_null_aa()
for src in ccat:
    ccat[src].compute(aa)
print "loading %d sources out of catalog (s) %s"%(len(ccat),','.join(ccats))


t = atpy.Table()
results = {}
for infile in args:
    print "working on: %s"%infile
    hdulist = pf.open(infile)
    header = hdulist[0].header
    proj = wcs.Projection(header)
    image = n.ma.masked_invalid(hdulist[0].data)
    print "WCS found types: ", proj.types
    print "WCS found units: ", proj.units
    try:
        map = proj.sub((proj.lonaxnum,proj.lataxnum))
    except:
        print "Aborting program. Could not find (valid) spatial map."
        raise
    for srcname,src in ccat.iteritems():
        #get the pixels within opts.radius of input source
        srcpos = n.array([src.ra,src.dec])*180/n.pi
        try:
            pix = query_disc(proj,srcpos,opts.radius)
            if pix.size<2:continue #skip if we don't get any pixels
        except:
            #print "source not in image"
            continue
        pmax = pix[image[pix[:,1],pix[:,0]].argmax(),:]
        try:
            RA,DEC = proj.toworld(pmax)
        except:
            continue
        print "looking for ",srcname,"@",srcpos
        Flux = image[pix[:,1],pix[:,0]].max()

        #get the pixels inside the surrounding annulus        
        doughpix = query_disc(proj,srcpos,outer_radius)
        holepix = query_disc(proj,srcpos,inner_radius)
        donutpix = n.array([px for px in doughpix if px not in holepix])
        RMS = n.std(image[donutpix[:,1],donutpix[:,0]])

        #record all the results
        results[srcname] = [RA,DEC,Flux,RMS]
        print n.array([RA,DEC,Flux,RMS])

#        cpix = proj.topixel(srcpos)
#        I = 20
#        print (cpix[0]-I,cpix[0]+I),(cpix[1]-I,cpix[1]+I)
#        ra = n.clip((cpix[0]-I,cpix[0]+I),0,proj.naxis[0])
#        dec = n.clip((cpix[1]-I,cpix[1]+I),0,proj.naxis[1])
#        print ra,dec
#        matshow(oimage[dec[0]:dec[1],
#                    ra[0]:ra[1]])
#        print image[dec[0]:dec[1],
#                    ra[0]:ra[1]].max()
#        show()
        #find the peak
        #return the coordinates of the peak
    outfile = '.'.join(infile.split('.')[:-1]+['vot'])
    names = results.keys()
    t.add_column('Name',names,dtype='S14',
            description='Name of source in finder catalog')
    t.add_column('RA',[results[name][0] for name in names],dtype='<f8',
            unit='deg',description='RA of peak')
    t.add_column('DEC',[results[name][1] for name in names],dtype='<f8',
            unit='deg',description='Dec of peak')
    t.add_column('Flux',[results[name][2] for name in names],dtype='<f8',
            unit='Jy',description='Magnitude of peak')
    t.add_column('RMS',[results[name][3] for name in names],dtype='<f8',
            unit='Jy',
            description='RMS in annulus between %d and %d of peak'%(inner_radius,outer_radius))
    print "Creating %s with %d sources"%(outfile,len(names))
    t.write(outfile)
