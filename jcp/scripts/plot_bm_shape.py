#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse
#from matplotlib.toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('bm_shape.py [options]')
a.scripting.add_standard_options(o, cal=True, src=True, pol=True, cmap=True)
o.add_option('-a', '--ant', dest='ant', type='int', default=0,
    help='Antenna to plot the beam for.  Default 0')
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='Frequency to plot beam for.  Default .150 GHz')
o.add_option('--res', dest='res', type='float', default=.01,
    help='Resolution of plot (in radians).  Default .01')
opts,args = o.parse_args(sys.argv[1:])

cmap = p.get_cmap(opts.cmap)
if not opts.src is None: srcs = opts.src.split(',')
else: srcs = []

aa = a.cal.get_aa(opts.cal, .1, opts.freq, 1)
print 'Plotting beam at %f GHz' % (opts.freq)
SZ = (int(n.pi/2/opts.res), int(2*n.pi/opts.res)+2)
lats,lons = n.indices(SZ)
lats = n.pi/2 - lats.astype(n.float) * opts.res
lons = lons.astype(n.float) * opts.res
#x,y,z = a.coord.azalt2top((lons.flatten(), lats.flatten()))
top = a.coord.azalt2top((lons.flatten(), lats.flatten()))
map = Basemap(projection='ortho', lat_0=90, lon_0=180)
map.drawmapboundary()
map.drawmeridians(n.arange(0, 360, 30))
map.drawparallels(n.arange(0, 90, 10))
cx, cy = map(lons * a.img.rad2deg, lats * a.img.rad2deg)
aa._cache = {'s_top':(0,0,1)}
aa.set_active_pol(opts.pol)
zresp = aa.bm_response(opts.ant, opts.ant)[:,0]
aa._cache = {'s_top':top}
data = aa.bm_response(opts.ant, opts.ant)[:,0] / zresp
#print data.shape, x.shape
#print data
data = n.clip(10*n.log10(n.abs(data)), -30, n.Inf)
data.shape = cx.shape
min = data.min()
max = data.max()
#max = data.max() / 2
step = (max - min) / 20
levels = n.arange(min-step, max+step, step)
map.contourf(cx, cy, data, levels, linewidths=0, cmap=cmap)
saz,salt = [], []
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
up = None
for t in n.arange(2454952.0, 2454953.0, .0001):
    aa.set_jultime(t)
    cat.compute(aa)
    _up = []
    for src in cat:
        src = cat[src]
        saz.append(src.az); salt.append(src.alt)
        if src.alt > .2: _up.append(src.src_name)
    _up.sort()
    if _up != up:
        up = _up
        print t, up
saz = n.array(saz); salt = n.array(salt)
#saz = n.concatenate([saz, saz+n.pi])
#salt = n.concatenate([salt, salt])
x,y = map(saz * a.img.rad2deg, salt * a.img.rad2deg)
map.plot(x, y, 'ko')
p.colorbar()
p.show()
