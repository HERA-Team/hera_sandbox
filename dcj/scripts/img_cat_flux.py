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
    #print len(n.logical_and(pix[:,0]<proj.naxis[0],pix[:,1]<proj.naxis[1]))
    return pix[n.logical_and(pix[:,0]<proj.naxis[0],pix[:,1]<proj.naxis[1])]
    #return S,R,D

clist,ccoff,ccats = a.scripting.parse_srcs(opts.src,opts.cat)
#print clist,ccoff,ccats
if not opts.cal is None:
    print "cal"
    ccat = a.cal.get_catalog(opts.cal,srcs=clist,cutoff=ccoff,catalogs=ccats)
else:
    ccat = a.src.get_catalog(srcs=clist,cutoff=ccoff,catalogs=ccats)
aa = a.scripting.get_null_aa()
for src in ccat:
    ccat[src].compute(aa)
print "loading %d sources out of catalog (s) %s"%(len(ccat),','.join(ccats))



for infile in args:
    t = atpy.Table()
    results = {}
    print "working on: %s"%infile
    hdulist = pf.open(infile)
    header = hdulist[0].header
    proj = wcs.Projection(header)
    image = n.ma.masked_invalid(hdulist[0].data).squeeze()
    iRA,iDec = n.mgrid[:image.shape[0],:image.shape[1]] #create an array to index the pixel locations
    print "WCS found types: ", proj.types
    print "WCS found units: ", proj.units
    freq = proj.sub((proj.specaxnum)).toworld(n.array([1]))[0] #grab the 
    try:
       MAP  = proj.sub((proj.lonaxnum,proj.lataxnum))
    except:
        print "Aborting program. Could not find (valid) spatial map."
        raise
    for srcname,src in ccat.iteritems():
        #get the pixels within opts.radius of input source
        P = src.get_params('*')
        ra,dec = n.degrees(P['ra']),n.degrees(P['dec'])
        catflux = P['jys']
        srcpos = n.array([ra,dec])
        try:
            pix = query_disc(MAP,srcpos,opts.radius)
            if pix.size<2 or pix[:,0].max()<0 or pix[:,1].max()<0:
                continue #skip if we don't get any pixels
        except:
            continue
        pmax = pix[image[pix[:,1],pix[:,0]].argmax(),:]
        #check that the pixel is on the image
        if n.any(pmax<0):
            continue
        try:
            RA,DEC = MAP.toworld(pmax)
        except:
            continue
        #find the peak in the subset
        Flux = image[pix[:,1],pix[:,0]].max()
        #find the location of the peak in the subset
        lmax = image[pix[:,1],pix[:,0]].argmax()
        lRA,lDec = iRA[pix[:,1],pix[:,0]].ravel()[lmax],iDec[pix[:,1],pix[:,0]].ravel()[lmax]
        pRA,pDec = MAP.toworld(n.array([lRA,lDec]))
        print "found ",srcname,"@",srcpos

        #get the pixels inside the surrounding annulus        
        doughpix = query_disc(MAP,srcpos,outer_radius)
        holepix = query_disc(MAP,srcpos,inner_radius)
        donutpix = n.array([px for px in doughpix if px not in holepix])
        RMS = n.std(image[donutpix[:,1],donutpix[:,0]])

        #compute the pixel number (raveled index for the facet image instead of healpix index) 
        srcpx = n.ravel_multi_index(pmax,dims=image.shape)
        #record all the results  
        #if n.any(n.isnan([Flux,RMS])):continue
        if not results.has_key(srcname): results[srcname] = []
        try:
            d = n.array([RA,DEC,pRA,pDec,RMS,srcpx,Flux,catflux,freq])
        except(ValueError):
            print [RA,DEC,pRA,pDec,RMS,srcpx,Flux,catflux,freq]
            raise
        TYPE = dtype([('RA','<f8',(1)),
               ('DEC','<f8'),
               ('pRA','<f8'),
               ('pDEC','<f8'),
               ('RMS','<f8'),
               ('srcpx','<i8'),
               ('Flux','<f8'),
               ('catflux','<f8'),
               ('freq','<f8')])
        results[srcname] = d.view(TYPE)
    if opts.outfile is None: outfile = '.'.join(infile.split('.')[:-1]+['vot'])
    else: outfile = opts.outfile
    names = results.keys()
    
    
    
    
    #save the data points in a vot file    
    t.add_column('srcpx',n.array([results[name]['srcpx'].squeeze() for name in names]),
            description='pixel number of peak.  do numpy.unravel_index(srcpx,dims=image.shape) to get pixel coords')
    print len(names)
    t.add_column('Name',names,dtype='S14',
            description='Name of source in %s catalog'%(opts.cat))
    t.add_column('RA',n.array([results[name]['RA'] for name in names]).squeeze(),dtype='<f4',
            unit='deg',description='RA of peak')
    t.add_column('DEC',n.array([results[name]['DEC'] for name in names]).squeeze(),dtype='<f4',
            unit='deg',description='Dec of peak')
    t.add_column('S_nu_',n.array([results[name]['Flux'] for name in names]).squeeze(),dtype='<f4',
            unit='Jy',description='Magnitude of peak')
    t.add_column('RMS',n.array([results[name]['RMS'] for name in names]).squeeze(),dtype='<f4',
            unit='Jy',
            description='RMS in annulus between %d and %d of peak'%(inner_radius,outer_radius))
    t.add_column('freq',n.array([results[name]['freq'] for name in names]).squeeze(),dtype='<f4',
            unit='MHz',
            description='Frequency in MHz of corresponding S_nu_ points.')        
    t.add_column('S_cat',n.array([results[name]['catflux'] for name in names]).squeeze(),dtype='<f4',
            unit='Jy',
            description='%s flux (extrapolated to 150MHz)'%(opts.cat))
    print "Creating %s with %d sources"%(outfile,len(names))
    t.write(outfile,overwrite=True)
