#!/usr/bin/env python
#
#  extract_subimg.py
#  
#
#  Created by Danny Jacobs on 1/19/10.
#  PAPER Project
#
"""
Given an input sources and fits image, extract a subimage.
WARNING/LIMITATIONS: Assumes square pixels, and doesnt work on cubes (yet)
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,ephem, pyfits as pf
from kapteyn import wcs
from pylab import *


o = optparse.OptionParser()
a.scripting.add_standard_options(o,src=True)
o.set_usage('extract_subimg.py -s srcs facets*fits')
o.add_option('--width',dest='width',default=10,type='float',
 help="Width of subimage in degrees [10]")
o.add_option('--sep',dest='sep',default=5,type='float',
   help='Selection criteria for source distance from img center. [5 degrees]')
o.add_option('-v',dest='verb',action='store_true',
    help='Print more')
o.add_option('--prefix',dest='prefix',
   help="Thumbnail filename prefix [optional]")
o.set_description(__doc__)
opts, args = o.parse_args(sys.argv[1:])
opts.sep *= a.img.deg2rad
def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)


srcs,coff,catalogs = a.scripting.parse_srcs(opts.src,opts.cat)
cat = a.src.get_catalog(srcs=srcs,catalogs=catalogs)
update_pos(cat)



for file in args:
    print file
    hdulist = pf.open(file)
    width = n.abs(opts.width/hdulist[0].header.get('CDELT1'))
    proj = wcs.Projection(hdulist[0].header)
    img_ra = hdulist[0].header.get('CRVAL1')*a.img.deg2rad
    img_dec = hdulist[0].header.get('CRVAL2')*a.img.deg2rad
    center = a.phs.RadioFixedBody(img_ra,
            img_dec)
    ephem.FixedBody.compute(center,ephem.J2000)
    for name,src in cat.iteritems():
        src_sep = ephem.separation(src,center)
        if src_sep<opts.sep:
            if opts.verb: print "getting \n",src
            ra = src.ra * a.img.rad2deg
            dec = src.dec * a.img.rad2deg
            if opts.verb:print src_sep,ra,dec
            px = proj.topixel((ra,dec,1,1))
            if opts.verb:print "at",px
            BL = n.round((px[0]-width/2,px[1]-width/2))
            TR = n.round((px[0]+width/2,px[1]+width/2))
            if opts.verb:print "subimg"
            if opts.verb: print BL, TR
            sub = n.meshgrid(range(BL[0],TR[0]), range(BL[1],TR[1]))
            im = n.flipud(hdulist[0].data.squeeze()[sub[1]+1,sub[0]+1])
#            imshow(im,aspect='equal')
            outfile = name+"_thumb.fits"
            if not opts.prefix is None: outfile = opts.prefix+outfile
            print ' > '+ outfile
            a.img.from_fits_to_fits(file,outfile,
                n.flipud(im),{'CRVAL1':ra,'CRVAL2':dec})