#!/usr/bin/env python
#
#  hpm_model.py
#  
#
#  Created by Danny Jacobs on 11/23/13.
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,healpy as hpy, ephem,os
import pickle
"""
Build a healpix map from a catalog.
example Usage:
hpm_model.py --nside=1024 -s 1/0.15 -C psa64_CSTbeam --cat=southern_sky_v3 -o mysky.fits
"""
def pix2radec(nside,pix):
    """input a pixel index into an nside healpix vector
        return ra,dec (in radians)
    """
    try:
        hp_loc = [hpy.pix2ang(nside,px) for px in pix]
    except(TypeError):
        hp_loc = [hpy.pix2ang(nside,pix)]
    ra = [loc[1] for loc in hp_loc]
    dec = [n.pi/2 - loc[0] for loc in hp_loc]
    return n.array(ra),n.array(dec)
def radec2pix(nside,ra,dec):
    if n.abs(dec)>n.pi or ra>2*n.pi: raise ValueError('Input coordinates must be in radians')
    theta = n.pi/2 - dec
    phi = ra
    return hpy.ang2px(nside,theta,phi)
def pixsep(nside,pix1,pix2):
    """return the seperation between two healpix pixels in degrees"""
    src1_loc = hpy.pix2ang(nside,pix1)
    src2_loc = hpy.pix2ang(nside, pix2)
    aa = a.scripting.get_null_aa()
    src1 = a.phs.RadioFixedBody(src1_loc[1],n.pi/2 - src1_loc[0])
    src2 = a.phs.RadioFixedBody(src2_loc[1],n.pi/2 - src2_loc[0])
    src1.compute(aa)
    src2.compute(aa)
    return ephem.separation(src1,src2)    
def pos2name(pos):
    raname=''.join(str(ephem.hours(pos['RA'])).split(':')).split('.')[0]
    decname=''.join(str(ephem.degrees(pos['DEC'])).split(':')).split('.')[0]
    if n.sign(pos['DEC'])>0.: decname = '+'+decname
    return raname+decname
def R(x):
    return n.round(x,2)
def Gaussian(x,width):
    return n.exp(-x**2/width**2)
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('-r',dest='radius',default=1.,type='float',
    help="psf calculation region around each input source, in degrees. [default = 1]")
o.add_option('--nside',default=1024,type=int,
    help="Nside of output map. default=1024")
o.add_option('-o',dest='output',default='catalog.fits',
    help='The output fits file.')
o.add_option('--osys',default='eq',
    help='Output coordinate system [eq(default),ga,ec]')
o.add_option('--psf',default=0.1,type=float,
    help='FWHM of source size when placed in map.[deg,default=0.1]')
#o.add_option('--hist',action='store_true',
#    help="Output the histograms of pixels in the rms annulus (pickled dict)")
opts, args = o.parse_args(sys.argv[1:])
outer_radius = opts.radius
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
n.set_printoptions(precision=2,suppress=True)
nside = opts.nside
sky = n.zeros((hpy.nside2npix(nside)))

for srcname,src in ccat.iteritems():
    if opts.osys=='ga':
        src = ephem.Galactic(src)
        catpx = hpy.ang2pix(nside,n.pi/2 -src.lat,src.lon)
        print src.lat,src.lon
    elif opts.osys=='ec':
        src = ephem.Equatorial(src)
        catpx = hpy.ang2pix(nside,n.pi/2 - src.lat,src.lon)
    else:
        catpx = hpy.ang2pix(nside,n.pi/2-src.dec,src.ra)
        
    v = hpy.pix2vec(nside,catpx)
    #account for different goddamn versions of healpy
    try: 
        outer_pix = hpy.query_disc(nside,
                    v,outer_radius,deg=True)
    except(TypeError):
        outer_pix = hpy.query_disc(nside,
                    v,outer_radius*n.pi/180)
        
    print "found %d pixels"%(len(outer_pix))
    #get the map pixels and calculate the psf
    print catpx, n.array(hpy.pix2ang(nside,catpx))*180/n.pi,ccat[srcname].jys.squeeze()
    _r = n.array([pixsep(nside,catpx,psfpix) for psfpix in outer_pix])
    psf = Gaussian(_r,opts.psf*n.pi/181)
    #insert the source
    flux = ccat[srcname].jys.squeeze()*psf
    if n.isnan(flux).sum()>0 or n.isnan(sky[outer_pix]).max():
        print 'NaN'
        continue
    sky[outer_pix] += flux
hpy.write_map(opts.output,sky,fits_IDL=False)

        
        
