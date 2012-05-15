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
o.add_option('--size', dest='size', type='str', default=60,
    help='Number of degrees on a side. [60_60]')
o.add_option('--res', dest='res', type='float', default=10,
    help='Size of pixel resolution in arcmin. [10]')
o.add_option('--isys', dest='isys', default='eq',
    help='Input coordinate system of (in map).  Can be eq (equatorial, default), ga (galactic), or ec (ecliptic).')
o.add_option('--osys', default='eq',
    help='Output coordinate system (in map).  Can be eq (equatorial, default) or ga (galactic).')
o.add_option('--fwidth', dest='fwidth', type='float', default=10,
    help='Width of gaussian, in degrees, used to downweight data away from pointing center for each facet.  Default 10.')
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data from healpix map.')
o.add_option('--ignore_weights',action='store_true',
    help='Dont use the weights column in the healpix fits table')
o.add_option('-j', '--juldate', dest='juldate', type='float', 
    help='Julian date used for locating moving sources.')
o.add_option('--output_weights',action='store_true',
    help='outputs a seperate weight map.')

opts, args = o.parse_args(sys.argv[1:])

def to_fits(filename,i,src,name,history=''):
#    filename = fname(ftag,cnt)
    print 'Saving data to', filename
    while len(i.shape) < 4: i.shape = i.shape + (1,)
#    cen = ephem.Equatorial(src._ra, src._dec)
#    # We precess the coordinates of the center of the image here to
#    # J2000, just to have a well-defined epoch for them.  For image coords to
#    # be accurately reconstructed, precession needs to be applied per pixel
#    # and not just per phase-center because ra/dec axes aren't necessarily
#    # aligned between epochs.  When reading these images, to be 100% accurate,
#    # one should precess the ra/dec coordinates back to the date of the
#    # observation, infer the coordinates of all the pixels, and then
#    # precess the coordinates for each pixel independently.
#    cen = ephem.Equatorial(cen, epoch=ephem.J2000)
    if opts.osys=='eq':
        axes = axes=('ra---sin', 'dec--sin')
        ra,dec = src._ra,src._dec
    elif opts.osys=='ga':
        axes = axes=('glon-sin','glat-sin')
        ra,dec = src.long,src.lat
    a.img.to_fits(filename, i, clobber=True,
        object=name, 
        ra=ra*a.img.rad2deg, dec=dec*a.img.rad2deg, epoch=2000.,
        d_ra=opts.res/60., d_dec=opts.res/60.,axes=axes,
        freq=150e6)

class Img:
    """
    A dumbed down version of the aipy Img class.  Impliments non-square dimensions
    but no gridding etc.
    Also the sizes are in actual image units, not uv. Minor difference.
    """
    def __init__(self, size=[100,100], res=1, mf_order=0):
        """size = the dimensions of the image in pixels
        res = The size of a (square) pixel in radians"""
        self.res = float(res)
        self.shape = size
    def get_LM(self, center=(0,0)):
        """Get the (l,m) image coordinates for an inverted UV matrix."""
        dim = self.shape
        M,L = n.indices(dim)
        L,M = n.where(L > dim[1]/2, -(L-dim[1]/2), (dim[1]/2 - L)), n.where(M > dim[0]/2, (M-dim[0]/2), -(dim[0]/2 - M))
        L,M = L.astype(n.float)*self.res, M.astype(n.float)*self.res
        mask = n.where(L**2 + M**2 >= 1, 1, 0)
        L,M = n.ma.array(L, mask=mask), n.ma.array(M, mask=mask)
        return L,M
    def get_top(self, center=(0,0)):
        """Return the topocentric coordinates of each pixel in the image."""
        x,y = self.get_LM(center)
        z = n.sqrt(1 - x**2 - y**2)
        return x,y,z
    def get_eq(self, ra=0, dec=0, center=(0,0)):
        """Return the equatorial coordinates of each pixel in the image, 
        assuming the image is centered on the provided ra, dec (in radians)."""
        x,y,z = self.get_top(center)
        shape,mask = x.shape, x.mask
        if len(mask.shape) == 0: mask = n.zeros(x.shape)
        vec = n.array([q.filled().flatten() for q in (x,y,z)])
        m = a.img.coord.top2eq_m(-ra, dec)
        vec = n.dot(m, vec)
        vec.shape = (3,) + shape
        return n.ma.array(vec, mask=[mask,mask,mask])    
if opts.osys=='ga':
    cat = {}
    for src in opts.src.split(','):
        cat[src] = ephem.Galactic(float(src.split('_')[0])*a.img.deg2rad,float(src.split('_')[1])*a.img.deg2rad)
elif not opts.cal is None:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
else:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.src.get_catalog(srclist, cutoff, catalogs)
o = ephem.Observer()
if opts.osys!='ga':
    if opts.juldate is None:
        o.date = ephem.J2000
        o.epoch = o.date
        try: del(cat['Sun'])
        except(KeyError): pass
    else:
        o.date = a.phs.juldate2ephem(opts.juldate)
        o.epoch = o.date
    for s in cat.values():
        try: a.phs.RadioFixedBody.compute(s, o)
        except(TypeError): a.phs.RadioSpecial.compute(s, o)

for infile in args:
    skymap = a.map.Map(fromfits=infile,interp=opts.interpolate)
    if opts.ignore_weights:
        print "ignoring weight column"
        Targs = ()
        Tkwargs = {'nside':skymap.nside(), 'scheme':skymap.scheme()}
        print Tkwargs
        skymap.wgt = a.healpix.HealpixMap(*Targs,**Tkwargs)
        skymap.wgt.set_map(n.ones_like(skymap.wgt.map))
    skymap.set_interpol(opts.interpolate)
    print "working on: ",infile
    for name,s in cat.iteritems():
        print "extracting: ",name
        if opts.osys=='ga':
            print "l = %s, b = %s"%(s.long,s.lat)
        else:
            print "RA = %s, DEC = %s"%(s._ra,s._dec)
        print opts.size,opts.res/60.
        RES = opts.res/60
        DIM = (n.array(opts.size.split('_')).astype(n.float)/RES)
        DIM = DIM.astype(n.int)
        RES *= a.img.deg2rad
        print "DIM (lat,long) = ",DIM,"  RES = ",RES
        im = Img(DIM, RES, mf_order=0)
        tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
        # Define a weighting for gridding data into the skymap
        map_wgts = n.exp(-(tx**2+ty**2) / n.sin(opts.fwidth*a.img.deg2rad)**2)
        map_wgts.shape = (map_wgts.size,)
        valid = n.logical_not(map_wgts.mask)
        map_wgts = map_wgts.compress(valid)
        if opts.osys=='ga':
            RA,DEC = s.long,s.lat
        else: RA,DEC = s._ra,s._dec
        x,y,z = im.get_eq(RA,DEC, center=(DIM/2,DIM/2))
        crd = n.row_stack((x.ravel(),y.ravel(),z.ravel()))
        print crd.shape
        m = a.coord.convert_m(opts.osys, opts.isys)
        print m.shape
        ex,ey,ez = n.dot(m, crd)
        ex = ex.compress(valid); ey = ey.compress(valid); ez = ez.compress(valid)
        img = skymap[ex,ey,ez]
        print img.shape
        img_wgts = skymap.wgt[ex,ey,ez]
        print img_wgts.shape
        print img.shape
        if opts.output_weights: 
            wgt_img = skymap.wgt[ex,ey,ez]
            wgt_img.shape = im.shape
            wgt_filename = infile[:-len('.fits')] +'_'+ name+'.wgt.fits'
            to_fits(wgt_filename,wgt_img,s,name,history='hpm_extract_facet:  Wgts for the tacet at %s extracted from healpix map %s [%s]'%(name,infile,time.asctime()))
        img.shape = im.shape
        filename = infile[:-len('.fits')] +'_'+ name+'.fits'
        to_fits(filename,img,s,name,history='hpm_extract_facet:  Facet at %s extracted from healpix map %s [%s]'%(name,infile,time.asctime()))
