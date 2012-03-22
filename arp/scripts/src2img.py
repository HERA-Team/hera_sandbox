#! /usr/bin/env python
import numpy as n, aipy as a, optparse, os, sys, ephem

o = optparse.OptionParser()
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('--center', dest='center', default='zen',
    help='Phase center of the image.  Default is zenith at the specified jultime')
o.add_option('--size', dest='size', type='float', default=400.,
    help='Size of maximum UV baseline.')
o.add_option('--res', dest='res', type='float', default=0.4,
    help='Resolution of UV matrix.')
o.add_option('--mfreq', dest='mfreq', default=.150, type='float',
    help='Frequency to generate image at (default is 0.150 GHz)')
o.add_option('-t', '--jultime', dest='jultime', type='float',
    help='Julian date to create zenith image at.')
o.add_option('-f', '--filename', dest='filename', 
    help='Filename to write to.')

opts, args = o.parse_args(sys.argv[1:])
aa = a.cal.get_aa(opts.cal, 1e-4, opts.mfreq, 1)

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
#a1s,a2s,ths = cat.get('srcshape')

aa.set_jultime(opts.jultime)
if opts.center == 'zen':
    cen = a.phs.RadioFixedBody(aa.sidereal_time(), aa.lat, epoch=aa.date, name='zen')
else:
    cen = a.cal.get_catalog(opts.cal, [opts.center], catalogs=catalogs).values()[0]
cat.compute(aa)
cen.compute(aa)

im = a.img.Img(size=opts.size, res=opts.res)
DIM = im.uv.shape[0]
imdata = n.zeros((DIM,DIM), dtype=n.float32)
ex,ey,ez = im.get_eq(ra=cen.ra, dec=cen.dec)
for src in cat:
    sx,sy,sz = cat[src].get_crds('eq', ncrd=3)
    i = n.argmin((ex-sx)**2 + (ey-sy)**2 + (ez-sz)**2)
    print src, i
    flx = cat[src].get_jys()
    imdata.flat[i] = flx

L,M = im.get_LM()
#for x,y,_sx,_sy,jys,name in zip(px,py,sx,sy,flx,cat.keys()):
#    print name, (x,y), (_sx,_sy), jys
#    print cat[name].ra, cat[name].dec
#    imdata[y,-x] = jys
#    #imdata[0,0] = jys
print n.where(a.img.recenter(imdata, (DIM/2,DIM/2)))

print 'Saving data to', opts.filename
cen = ephem.Equatorial(cen, epoch=ephem.J2000)
a.img.to_fits(opts.filename, a.img.recenter(imdata, (DIM/2,DIM/2)), clobber=True,
    object=opts.center, obs_date=str(aa.date),
    ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
    d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[1,1]*a.img.rad2deg,
        freq=opts.mfreq)

