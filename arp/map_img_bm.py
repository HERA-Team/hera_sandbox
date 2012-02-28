#! /usr/bin/env python
"""
This is a general-purpose script for making images from MIRIAD UV files.  Data
(optionally selected for baseline, channel) are read from the file, phased
to a provided position, normalized for passband/primary beam effects, gridded
to a UV matrix, and imaged
"""

import aipy as a, numpy as n, pylab as p
import sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('map_img_bm.py [options] img1.fits ...')
o.set_description(__doc__)
#a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('-m', '--map', dest='map', 
    help='Map FITS filename.')
opts, args = o.parse_args(sys.argv[1:])


h = a.map.Map(fromfits=opts.map)
h.set_interpol(True)

prev_dra = None
bmsum,bmwgt = 0,0
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
    if prev_dra != kwds['d_ra']:
        prev_dra = kwds['d_ra']
        DIM = img.shape[0]
        RES = 1. / (kwds['d_ra'] * a.img.deg2rad * DIM)
        im = a.img.Img(DIM*RES, RES, mf_order=0)
        tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
        L,M = im.get_LM()
        mask = tx.mask
    # Get coordinates of image pixels in the epoch of the observation
    ex,ey,ez = im.get_eq(ra, dec, center=(DIM/2,DIM/2))
    ex = ex.filled(0).flatten(); ey = ey.filled(0).flatten(); ez = ez.filled(0).flatten()
    # Precess the pixel coordinates to the (J2000) epoch of the map
    m = a.coord.convert_m('eq','eq', iepoch=kwds['obs_date'], oepoch=ephem.J2000)
    ex,ey,ez = n.dot(m, n.array([ex,ey,ez]))
    d = h[ex,ey,ez]
    d.shape = mask.shape
    d = n.where(mask, 0, d)
    if True:
        w = h.wgt[ex,ey,ez] # weighting factor relating to trustworthiness of map
        w.shape = mask.shape
        w = n.where(mask, 0, w)
        w = w.clip(0,w.max()/1e2)
    else: w = 1
    bmsum += img * w**2 *n.abs(d)
    bmwgt += n.abs(d)**2 * w**2
    if False:
        if i == 0:
            p.ion()
            p.subplot(131); plt1 = p.imshow(img, vmax=10, vmin=-1, origin='lower', interpolation='nearest')
            p.subplot(132); plt2 = p.imshow(  d, vmax=10, vmin=-1, origin='lower', interpolation='nearest')
            p.subplot(133); plt3 = p.imshow(bmsum/n.where(bmwgt > 0, bmwgt, 1e6), vmax=2, vmin=0, origin='lower', interpolation='nearest')
        else:
            plt1.set_data(img)
            plt2.set_data(d)
            plt3.set_data(bmsum / n.where(bmwgt > 0, bmwgt, 1e6))
        p.draw()

def to_fits(i,filename):
    print 'Saving data to', filename
    while len(i.shape) < 4: i.shape = i.shape + (1,)
    a.img.to_fits(filename, i, clobber=True,
        object='zen', obs_date=str(ephem.now()),
        ra=0, dec=90, epoch=2000.,
        d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[1,1]*a.img.rad2deg,
        freq=.150,history='')

to_fits(bmsum, 'beam.dim.fits')
to_fits(bmwgt, 'beam.dbm.fits')

