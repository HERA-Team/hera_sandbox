#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse
from mpl_toolkits.basemap import Basemap

p.ion()

o = optparse.OptionParser()
o.set_usage('skyview.py [options]')
a.scripting.add_standard_options(o, cal=True, src=True, cmap=True, max=True, drng=True)
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='Frequency to plot beam for.  Default .150 GHz')
o.add_option('-j', '--juldate', dest='juldate', type='float',
    help='Julian date used for locating moving sources.')
o.add_option('--src_mark', dest='src_mark', default='',
    help='Marker to put on src locations.  Can be: ".,o,+,x,^,v".  Default no marker.')
o.add_option('--src_color', dest='src_color', default='k',
    help='Color of source label.  Can be: "k,w,r,b".  Default "k".')
o.add_option('-o', '--outfile', dest='outfile', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')
o.add_option('--isys', dest='isys', default='eq',
    help='Input coordinate system (in map).  Can be eq (equatorial, default), ga (galactic), or ec (ecliptic).')
o.add_option('--iepoch', dest='iepoch', type='float', default=a.ephem.J2000,
    help='Epoch of input coordinates (in map).  Default J2000.')
o.add_option('--nobar', dest='nobar', action='store_true',
    help="Do not show colorbar.")
o.add_option('--res', dest='res', type='float', default=0.25,
    help="Resolution of plot (in degrees).  Default 0.25.")
o.add_option('--nside', dest='nside', type='int',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
opts,args = o.parse_args(sys.argv[1:])

cmap = p.get_cmap(opts.cmap)
if not opts.src is None: srcs = opts.src.split(',')
else: srcs = []

aa = a.cal.get_aa(opts.cal, .1, opts.freq, 1)
aa.set_jultime(opts.juldate)
print 'Plotting sky at %f GHz' % (opts.freq)

# Create plotting grid
#SZ = (int(n.pi/2/opts.res), int(2*n.pi/opts.res)+2)
#lats,lons = n.indices(SZ)
#lats = n.pi/2 - lats.astype(n.float) * opts.res
#lons = lons.astype(n.float) * opts.res
#top = a.coord.azalt2top((lons.flatten(), lats.flatten()))
#map = Basemap(projection='ortho', lat_0=90, lon_0=180)
map = Basemap(projection='ortho', lat_0=aa.lat, lon_0=-a.img.rad2deg*aa.sidereal_time())
map.drawmapboundary()
map.drawmeridians(n.arange(0, 360, 30))
map.drawparallels(n.arange(-90, 90, 30))
#cx, cy = map(lons * a.img.rad2deg, lats * a.img.rad2deg)
#aa._cache = {'s_top':(0,0,1)}
#zresp = aa.bm_response(0, 0, pol='yy')[:,0]
#aa._cache = {'s_top':top}
#bm = aa.bm_response(0, 0, pol='yy')[:,0] / zresp
#bm = n.clip(10*n.log10(n.abs(bm)), -100, n.Inf)
#bm.shape = cx.shape
#min = -30.
#max = 0.
#step = (max - min) / 10
#levels = n.arange(min-step, max+step, step)

print 'Reading %s' % args[0]
h = a.map.Map(fromfits=args[0])
print 'SCHEME:', h.scheme()
print 'NSIDE:', h.nside()
if not opts.nside is None:
    nh = a.healpix.HealpixMap(nside=opts.nside)
    nh.from_hpm(h)
    h = nh
h.set_interpol(False)

# Convert plot coords into equatorial
img = a.img.Img(300,.5)
ex,ey,ez = img.get_eq(aa.sidereal_time(), aa.lat)
tx,ty,tz = img.get_top()
shape = ex.shape
ex,ey,ez = ex.flatten(), ey.flatten(), ez.flatten()
tx,ty,tz = tx.flatten(), ty.flatten(), tz.flatten()
aa._cache = {'s_top':(0,0,1)}
zresp = aa.bm_response(0, 0, pol='yy')[:,0]
aa._cache = {'s_top':(tx,ty,tz)}
bm = aa.bm_response(0, 0, pol='yy')[:,0] / zresp
bm.shape = shape
bm = a.img.recenter(bm, n.array(bm.shape)/2)

try: data, indices = h[ex,ey,ez]
except(ValueError): data = h[ex,ey,ez]
data.shape = shape
data = a.img.recenter(data, n.array(data.shape)/2)

data = n.log10(n.abs(bm*data))
max = opts.max
if not max: max = data.max()
if not opts.drng: min = data.min()
else: min = max - opts.drng
data = data.clip(min,max)

plt = {}
if opts.src:
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
    cat.compute(aa)
    for src in cat.values():
        if src.alt > 0:
            x,y = map(src.az * a.img.rad2deg, src.alt * a.img.rad2deg)
        else: x,y = -1,-1
        plt[src.src_name] = p.plot(x,y, 'ko')[-1]

plt['im'] = map.imshow(data, vmax=max, vmin=min, cmap=cmap, interpolation='nearest')
#map.contour(cx, cy, bm, levels, cmap=cmap)
if not opts.nobar: p.colorbar(shrink=.5)

if False:
  for jd in n.arange(opts.juldate, opts.juldate+1, .01):
    print jd
    aa.set_jultime(jd)
    ex,ey,ez = img.get_eq(aa.sidereal_time(), aa.lat)
    ex,ey,ez = ex.flatten(), ey.flatten(), ez.flatten()
    try: data, indices = h[ex,ey,ez]
    except(ValueError): data = h[ex,ey,ez]
    data.shape = shape
    data = a.img.recenter(data, n.array(data.shape)/2)
    data = n.clip(n.log10(n.abs(bm*data)), min, max)
    plt['im'].set_data(data)
    if opts.src:
        cat.compute(aa)
        for src in cat.values():
            if src.alt > 0:
                x,y = map(src.az * a.img.rad2deg, src.alt * a.img.rad2deg)
                plt[src.src_name].set_xdata(x)
                plt[src.src_name].set_ydata(y)
    p.draw()

p.show()
