#!/usr/bin/env python
#
#  sky_gen.py
#  
#
#  Created by Danny Jacobs on 6/3/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
import aipy as a, numpy as n, pylab as p, sys, os, ephem, optparse

o = optparse.OptionParser()
o.set_usage('sky_gen.py [options] ')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True)
o.add_option('--nside', dest='nside', type='int',default=1024,
    help="Set NSIDE to a power of 2.  Default=1024")
o.add_option('-o', '--outfile', dest='outfile', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')
o.add_option('-w','--width',dest='scale',default=20,type='float',
   help="""The size of a source in arcmin. Default=20.""")
o.add_option('-r','--res',dest='res',default=2,type='float',
   help="""The resolution of a source in arcmin. Default=2.""")
opts,args = o.parse_args(sys.argv[1:])


def _gen_K(sigma):
    """
    Generate a gaussian kernel matrix of FWHM=sigma [columns]
    """
    dim = int(sigma*12.0)
    x = n.linspace(-dim/2,dim/2,num=dim)
    X,Y = n.meshgrid(x,x)
    R = n.sqrt(X**2+Y**2)
    return dim, 1/n.sqrt(2*n.pi*sigma**2)*n.exp(-(R**2.0)/(2*sigma**2))

dim, K = _gen_K(opts.scale/opts.res)
skymap = a.map.Map(nside=opts.nside)
skymap.set_interpol(True)
size = dim*opts.res*a.const.arcmin
print "Adding sources with dim = ",dim
if not opts.src is None:
    srclist,cutoff = a.scripting.parse_srcs(opts.src)
    cat = a.src.get_catalog(srcs=srclist, cutoff=cutoff)
    o = ephem.Observer()
    for s in cat.values():
        try: a.phs.RadioFixedBody.compute(s, o)
        except(TypeError): a.phs.RadioSpecial.compute(s, o)
    scrds = [ephem.Equatorial(s.ra,s.dec,epoch=ephem.J2000) for s in cat.values()]
    print "Adding ",str(len(scrds))," sources"
    print "Generating ",str(len(scrds)*dim*dim)," pixels"
    afreqs = n.array([.150])
    cat.update_jys(afreqs)
    sflxs = cat.get_jys().squeeze()
    try: sflxs[0]
    except(IndexError):sflxs=[sflxs]
    sflxs = n.clip(sflxs,10,10e10)
    snams = cat.keys()
    for i,s in enumerate(scrds):
        ra = n.linspace(s.ra-size/2.0,s.ra+size/2.0,num=dim)
        dec = n.linspace(s.dec-size/2.0,s.dec+size/2.0,num=dim)
        RA,DEC = n.meshgrid(ra,dec)
#        print n.max(sflxs[i]*K.flatten())
        skymap.add((DEC.flatten()-n.pi/2,RA.flatten()),n.ones_like(RA.flatten()),sflxs[i]*K.flatten())
skymap.to_fits(opts.outfile,clobber=True)
#a.img.to_fits('src_img1.fits',K,axes=('ra---sin','dec---sin'),d_ra=0.1,d_dec=0.1,ra=180)
#mk_map.py src_img1.fits -m test_src5.fits
#skymap = a.map.Map(nside=opts.nside)
#skymap.set_interpol(opts.interpolate)
#
## Plot src labels and markers on top of map image
#if not opts.src is None:
#    sx, sy = map(slons,slats)
#    for name, xpt, ypt, flx in zip(snams, sx, sy, sflxs):
#        if xpt >= 1e30 or ypt >= 1e30: continue
#        if opts.src_mark != '':
#            map.plot(sx, sy, opts.src_color+opts.src_mark,markerfacecolor=None)
#        if flx < 10: flx = 10
#        p.text(xpt+.001, ypt+.001, name, size=5+2*int(n.round(n.log10(flx))),
#            color=opts.src_color)
#
#print 'Saving to', opts.outfile
#p.savefig(opts.outfile)