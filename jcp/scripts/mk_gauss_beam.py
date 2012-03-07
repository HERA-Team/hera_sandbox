#!/usr/global/paper/bin/python

import aipy as a, numpy as n, pylab as p, sys, os, optparse
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_beam.py [options] mapfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cmap=True, max=True, drng=True,cal=True)
o.add_option('--npix', dest='npix',type='int', default=1000,
    help="Number of pixels used in the sky image.")
o.add_option('--sigx',dest='sigx',type='float',default=.5,
    help="Standard deviation of the Guassian in the x direction.")
o.add_option('--sigy',dest='sigy',type='float',default=.5,
    help="Standard deviation of the Guassian in the y direction.")

opts,args = o.parse_args(sys.argv[1:])
center = opts.npix/2

im = a.img.Img(size=(opts.npix*.4),res=.4)
x,y,z = im.get_top()

beam = n.exp(-(x/opts.sigx)**2)*n.exp(-(y/opts.sigy)**2)
mbeam = n.ma.array(beam,mask=x.mask)

outname = 'gbeam%sx%sy' % (opts.sigx, opts.sigy) + '.fits'
outfile = a.map.Map(nside=256)
xx,yy,zz = x.compressed(),y.compressed(),z.compressed()
hbeam = mbeam.compressed()
hbeam.shape = xx.shape
outfile.add((xx,yy,zz),1,hbeam)
outfile.to_fits(outname,clobber=True)
