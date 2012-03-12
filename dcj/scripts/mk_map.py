#! /usr/bin/env python
"""
This is a general-purpose script for making faceted, spherical maps (stored in
Healpix FITS format) from individual "flat" maps stored in FITS files.

Author: Aaron Parsons
"""

import sys, numpy as n, os, aipy as a, optparse, ephem,pyfits as pf

o = optparse.OptionParser()
o.set_usage('mk_map.py [options] *.fits')
o.set_description(__doc__)
o.add_option('-i', '--interpolate', dest='interpolate', action='store_true',
    help='Use sub-pixel interpolation when gridding data to healpix map.')
o.add_option('-m', '--map', dest='map',
    help='The skymap file to use.  If it exists, new data will be added to the map.  Othewise, the file will be created.')
o.add_option('--nside', dest='nside', type='int', default=256,
    help='NSIDE parameter for map, if creating a new file.')
o.add_option('--fwidth', dest='fwidth', type='float', default=10,
    help='Width of gaussian, in degrees, used to downweight data away from pointing center for each facet.  Default 10.')
o.add_option('--crop',type='float',default=None,
    help='Optional radius by which to crop input images.')
#TODO add crop option.
opts, args = o.parse_args(sys.argv[1:])
if not opts.crop is None:
    print "cropping input images at a radius of %5.2f"%opts.crop
# Open skymap
if os.path.exists(opts.map): skymap = a.map.Map(fromfits=opts.map)
else: skymap = a.map.Map(nside=opts.nside)
skymap.set_interpol(opts.interpolate)

prev_dra = None
for i, filename in enumerate(args):
    img, kwds = a.img.from_fits(filename)
    hdulist = pf.open(filename)
    history = hdulist[0].header.get_history()
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
    if prev_dra != kwds['d_ra']:
        prev_dra = kwds['d_ra']
        DIM = img.shape[0]
        RES = 1. / n.abs(kwds['d_ra'] * a.img.deg2rad * DIM)
        im = a.img.Img(DIM*RES, RES, mf_order=0)
        tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
        # Define a weighting for gridding data into the skymap
        map_wgts = n.exp(-(tx**2+ty**2) / n.sin(opts.fwidth*a.img.deg2rad)**2)
        if not opts.crop is None:
            R = n.sqrt(tx**2 + ty**2)
            crop_mask = n.where(R>n.sin(opts.crop*a.img.deg2rad),1,0)
            map_wgts.mask = n.logical_or(map_wgts.mask,crop_mask.data)
        map_wgts.shape = (map_wgts.size,)
        valid = n.logical_not(map_wgts.mask)
        map_wgts = map_wgts.compress(valid)
    # Get coordinates of image pixels in the epoch of the observation
    ex,ey,ez = im.get_eq(ra, dec, center=(DIM/2,DIM/2))
    ex = ex.compress(valid); ey = ey.compress(valid); ez = ez.compress(valid)
    img = img.flatten(); img = img.compress(valid)
    # Precess the pixel coordinates to the (J2000) epoch of the map
    m = a.coord.convert_m('eq','eq', 
        iepoch=kwds['obs_date'], oepoch=ephem.J2000)
    ex,ey,ez = n.dot(m, n.array([ex,ey,ez])) 
    # Put the data into the skymap
    skymap.add((ex,ey,ez), map_wgts, img)
if len(history)==0: history = ['Warning: No facet history.']
history =  '\n'.join([h.strip() for h in history])+ '\n' +\
            ' '.join(sys.argv) +'\n'
print history
skymap.to_fits(opts.map, clobber=True,history=history)

