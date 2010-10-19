#!/usr/bin/env python
#
#  hpm_extract_facet.py
#  
#
#  Created by Danny Jacobs on 9/7/10.
#  PAPER Project
#
"""
Extracts a facet from a healpix map.  Makes a small angle approximation,
so don't try this with images larger than ~15 degrees.
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,time,ephem

o = optparse.OptionParser()
a.scripting.add_standard_options(o, src=True,cal=True)
o.set_description(__doc__)
o.add_option('--size', dest='size', type='float', default=60,
    help='Number of degrees on a side. [60]')
o.add_option('--res', dest='res', type='float', default=10,
    help='Size of pixel resolution in arcmin. [10]')
o.add_option('--isys', dest='isys', default='eq',
    help='Input coordinate system (in map).  Can be eq (equatorial, default), ga (galactic), or ec (ecliptic).')
o.add_option('--fwidth', dest='fwidth', type='float', default=10,
    help='Width of gaussian, in degrees, used to downweight data away from pointing center for each facet.  Default 10.')
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data from healpix map.')
opts, args = o.parse_args(sys.argv[1:])

def to_fits(filename,i,src,history=''):
#    filename = fname(ftag,cnt)
    print 'Saving data to', filename
    while len(i.shape) < 4: i.shape = i.shape + (1,)
    cen = ephem.Equatorial(src._ra, src._dec)
    # We precess the coordinates of the center of the image here to
    # J2000, just to have a well-defined epoch for them.  For image coords to
    # be accurately reconstructed, precession needs to be applied per pixel
    # and not just per phase-center because ra/dec axes aren't necessarily
    # aligned between epochs.  When reading these images, to be 100% accurate,
    # one should precess the ra/dec coordinates back to the date of the
    # observation, infer the coordinates of all the pixels, and then
    # precess the coordinates for each pixel independently.
    cen = ephem.Equatorial(cen, epoch=ephem.J2000)
    a.img.to_fits(filename, i, clobber=True,
        object=src.src_name, 
        ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
        d_ra=opts.res/60., d_dec=opts.res/60.,
        freq=150e6,history=history)


if not opts.cat is None:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
else:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.src.get_catalog(opts.cal, srclist, cutoff, catalogs)


#infile = args[0]
for infile in args:
    skymap = a.map.Map(fromfits=infile,interp=opts.interpolate)
    skymap.set_interpol(opts.interpolate)
    print "working on: ",infile
    for s in cat.values():
        print "extracting: ",s.src_name
        print "RA = %s, DEC = %s"%(s._ra,s._dec)
        print opts.size,opts.res/60.
        DIM = int(opts.size/(opts.res/60.))
        RES = int(1. / (opts.size * a.img.deg2rad))
        print "DIM = ",DIM,"  RES = ",RES
        im = a.img.Img(DIM*RES, RES, mf_order=0)
        tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
        print tx
        # Define a weighting for gridding data into the skymap
        map_wgts = n.exp(-(tx**2+ty**2) / n.sin(opts.fwidth*a.img.deg2rad)**2)
        map_wgts.shape = (map_wgts.size,)
        valid = n.logical_not(map_wgts.mask)
        map_wgts = map_wgts.compress(valid)
        if opts.isys=='ga': 
            RA,DEC = ephem.Galactic(ephem.Equatorial(s._ra,s._dec)).long,\
                    ephem.Galactic(ephem.Equatorial(s._ra,s._dec)).lat
            print RA,DEC
        else: RA,DEC = s._ra,s._dec
        ex,ey,ez = im.get_eq(RA,DEC, center=(DIM/2,DIM/2))
        ex = ex.compress(valid); ey = ey.compress(valid); ez = ez.compress(valid)
        print ex
        img = skymap[ex,ey,ez]
        img.shape = im.shape
        filename = infile[:-len('.fits')] +'_'+ s.src_name+'.fits'
        to_fits(filename,img,s,history='hpm_extract_facet:  Facet at %s extracted from healpix map %s [%s]'%(s.src_name,infile,time.asctime()))