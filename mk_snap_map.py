#!/usr/global/paper/bin/python
"""
This is a general-purpose script for making faceted, spherical maps (stored in
Healpix FITS format) from individual "flat" maps stored in FITS files.

Author: Aaron Parsons
"""

import sys, numpy as n, os, aipy as a, optparse, ephem

o = optparse.OptionParser()
o.set_usage('mk_map.py [options] *.fits')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data to healpix map.')
o.add_option('-m', '--map', dest='map',
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created.')
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='The frequency to use for beam weighting')
o.add_option('--nside', dest='nside', type='int', default=256,
    help='NSIDE parameter for map, if creating a new file.')
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, .001, opts.freq, 1)

# Open skymap
if os.path.exists(opts.map): skymap = a.map.Map(fromfits=opts.map)
else: skymap = a.map.Map(nside=opts.nside)
skymap.set_interpol(opts.interpolate)

bm_wgts = None
for i, filename in enumerate(args):
    img, kwds = a.img.from_fits(filename)
    img = img.squeeze()
    # Read ra/dec of image center, which are stored in J2000
    assert(kwds['epoch'] == 2000)
    s = ephem.Equatorial(kwds['ra']*a.img.deg2rad, kwds['dec']*a.img.deg2rad, 
        epoch=ephem.J2000)
    # To precess the entire image to J2000, we actually need to precess the
    # coords of the center back to the epoch of the observation (obs_date),
    # get the pixel coordinates in that epoch, and then precess the coords of
    # each pixel.  This is because extrapolating pixel coords from the J2000 
    # center assumes the J2000 ra/dec axes, which may be tilted relative to
    # the ra/dec axes of the epoch of the image.
    s = ephem.Equatorial(s, epoch=kwds['obs_date'])
    ra, dec = s.get()
    print '-----------------------------------------------------------'
    print 'Reading file %s (%d / %d)' % (filename, i + 1, len(args))
    print kwds
    print 'Pointing (ra, dec):', ra, dec
    print 'Image Power:', n.abs(img).sum()
    if bm_wgts is None:
        print 'Calculating beam weighting'
        DIM = img.shape[0]
        RES = 1. / n.abs(kwds['d_ra'] * a.img.deg2rad * DIM)
        im = a.img.Img(DIM*RES, RES, mf_order=0)
        tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
        valid = n.logical_not(tx.mask).flatten()
        tx = tx.flatten().compress(valid)
        ty = ty.flatten().compress(valid)
        tz = tz.flatten().compress(valid)
        bm_wgts = aa[0].bm_response((tx,ty,tz), pol=opts.pol[0])
        if opts.pol[0] == opts.pol[1]: bm_wgts *= bm_wgts
        else: bm_wgts *= aa[0].bm_response((tx,ty,tz), pol=opts.pol[1])
        bm_wgts = n.abs(bm_wgts.squeeze())
        bm_wgts_clip = n.where(bm_wgts == 0, 1, bm_wgts)
        map_wgts = bm_wgts**2
        print 'Done calculating beam weighting'
    # Get coordinates of image pixels in the epoch of the observation
    ex,ey,ez = im.get_eq(ra, dec, center=(DIM/2,DIM/2))
    ex = ex.flatten().compress(valid)
    ey = ey.flatten().compress(valid)
    ez = ez.flatten().compress(valid)
    # Precess the pixel coordinates to the (J2000) epoch of the map
    print 'Precessing coords'
    m = a.coord.convert_m('eq','eq', 
        iepoch=kwds['obs_date'], oepoch=ephem.J2000)
    ex,ey,ez = n.dot(m, n.array([ex,ey,ez])) 
    img = img.flatten().compress(valid) / bm_wgts_clip  # Remove beam response
    if False: # Try to guess flux scale of image (doesn't work)
        print 'Guessing flux scale'
        flx_wgts = map_wgts * skymap.wgt[ex,ey,ez]
        flx_wgts_sum = flx_wgts.sum()
        if flx_wgts_sum == 0: gain = 1.
        else:
            skymap_flx = n.abs(skymap[ex,ey,ez]) * flx_wgts
            img_flx = n.abs(img) * flx_wgts
            gain = img_flx.sum() / skymap_flx.sum()
            print 'Gain:', gain, img_flx.sum(), skymap_flx.sum()
        print 'Gain:', gain
        img /= gain
    # Put the data into the skymap
    skymap.add((ex,ey,ez), map_wgts, img) # down-weight by beam resp squared
skymap.to_fits(opts.map, clobber=True)

