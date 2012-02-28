#! /usr/bin/env python
"""
This is a general-purpose script for making images from MIRIAD UV files.  Data
(optionally selected for baseline, channel) are read from the file, phased
to a provided position, normalized for passband/primary beam effects, gridded
to a UV matrix, and imaged
"""

import aipy as a, numpy as n, sys, optparse, ephem, os

o = optparse.OptionParser()
o.set_usage('map2img.py [options] mapfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('-o', '--output', dest='outfile', default='out.fits',
    help='Output filename.  Default out.fits')
o.add_option('--size', dest='size', type='int', default=300,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--beam_pol', dest='beam_pol', 
    help='Polarization (xx,yy,xy,yx) of the primary beam to apply to the image (i.e. make it a perceived Jy image)')
o.add_option('--mfreq', dest='mfreq', type='float', default=.150,
    help='Create an image evaluated at this frequency.  Default is .150 GHz')
opts, args = o.parse_args(sys.argv[1:])

def to_fits(i,src,history=''):
    print 'Saving data to', opts.outfile
    while len(i.shape) < 4: i.shape = i.shape + (1,)
    cen = ephem.Equatorial(src.ra, src.dec, epoch=aa.epoch)
    # We precess the coordinates of the center of the image here to
    # J2000, just to have a well-defined epoch for them.  For image coords to
    # be accurately reconstructed, precession needs to be applied per pixel
    # and not just per phase-center because ra/dec axes aren't necessarily
    # aligned between epochs.  When reading these images, to be 100% accurate,
    # one should precess the ra/dec coordinates back to the date of the
    # observation, infer the coordinates of all the pixels, and then
    # precess the coordinates for each pixel independently.
    cen = ephem.Equatorial(cen, epoch=ephem.J2000)
    a.img.to_fits(opts.outfile, i, clobber=True,
        object=src.src_name, obs_date=str(aa.date),
        ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
        d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[1,1]*a.img.rad2deg,
        freq=n.average(aa[0].beam.afreqs),history=history)

aa = a.cal.get_aa(opts.cal, n.array([opts.mfreq]))

# Get all sources that will be used as phase centers. 
if opts.src == 'zen':
    srcs = [a.phs.RadioFixedBody(aa.sidereal_time(), aa.lat, name='zen')]
    cat = a.phs.SrcCatalog(srcs)
elif not opts.src is None: 
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)

h = a.map.Map(fromfits=args[0])
h.set_interpol(True)

DIM = int(opts.size/opts.res)
for s in cat:
    src = cat[s]
    src.compute(aa)
    im = a.img.Img(opts.size, opts.res, mf_order=0)
    L,M = im.get_LM()
    if opts.beam_pol != None:
        p1,p2 = opts.beam_pol
        tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
        SH = tx.shape
        tx,ty,tz = tx.flatten(), ty.flatten(), tz.flatten()
        bm_resp = aa[0].bm_response((tx,ty,tz), pol=p1) * n.conj(aa[0].bm_response((tx,ty,tz), pol=p2))
        bm_resp.shape = SH
    else: bm_resp = 1
    x,y,z = im.get_eq(ra=src._ra, dec=src._dec, center=(DIM/2,DIM/2))
    mask = x.mask
    d = h[x.filled(0).flatten(), y.filled(0).flatten(), z.filled(0).flatten()]
    d.shape = x.shape
    d = n.where(mask, 0, d * bm_resp)
    to_fits(d,src)

