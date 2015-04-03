#!/usr/global/paper/bin/python

import aipy as a, numpy as n, pylab as p, sys, optparse, ephem

o = optparse.OptionParser()
o.add_option('--size', dest='size', type='int', default=400,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.5,
    help='Resolution of UV matrix.')
o.add_option('--no_w', dest='no_w', action='store_true',
    help="Don't use W projection.")
o.add_option('--wres', dest='wres', type='float', default=0.5,
    help="W-Plane projection resolution.  Default 0.5")
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

binsize = 4
nbins = n.float(len(args)/binsize)

ind = 1
rmap,tmap = 0,0
for map in args:
    map = a.img.from_fits(map)[0]
    if ind < binsize:
        #print 'in bin', ind 
        tmap += map
        ind += 1
    elif ind == binsize:
        #print 'bin done', ind
        rmap += tmap**2
        tmap = 0
        ind = 1

rmap = (rmap/nbins)**.5

#p.imshow(rmap,interpolation='nearest')
#p.show()

#fits stuff
filename = 'rmap.fits'
aa = a.cal.get_aa(opts.cal, n.array([.150]))
ra,dec = '21:52:10.67','-30:43:17.5'
src = a.phs.RadioFixedBody(ra,dec)
src.compute(aa)
cen = ephem.Equatorial(src.ra, src.dec, epoch=aa.epoch)
cen = ephem.Equatorial(cen, epoch=ephem.J2000)

if opts.no_w:
    im = a.img.Img(opts.size, opts.res, mf_order=0)
else:
    im = a.img.ImgW(opts.size, opts.res, mf_order=0, wres=opts.wres)
L,M = im.get_LM()


a.img.to_fits(filename, rmap, clobber=True,
    #object=src.src_name, obs_date=str(aa.date),
    ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
    d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[1,1]*a.img.rad2deg,
    freq=n.average(aa[0].beam.afreqs),history='')

