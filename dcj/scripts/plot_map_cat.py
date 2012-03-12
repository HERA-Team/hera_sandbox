#!/usr/global/paper/bin/python
"""
Script for displaying a projection of a spherical (Healpix) data set stored
in a *.fits file.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, pylab as p, sys, os, ephem, optparse

class Basemap:
    """A placeholder class to give plot_map.py some functionality if
    matplotlib-basemap is not installed."""
    def __init__(self, projection, lat_0=0, lon_0=0, **kwargs):
        if projection != 'cyl':
            raise ValueError('Without matplotlib-basemap installed, only cyl projection is supported.')
        self.lat_0, self.lon_0 = lat_0, lon_0
    def __call__(self, lon, lat, inverse=False):
        if inverse:
            lon = lon.astype(n.float) * 90 - self.lon_0
            lon = self.wrap(lon, -180, 180)
            lat = lat.astype(n.float) * 90 - self.lat_0
            lat = self.wrap(lat, -90, 90)
            return lon,lat
        else:
            x = (lon + self.lon_0) / 90.
            x = self.wrap(x, -2, 2)
            y = (lat + self.lat_0) / 90.
            y = self.wrap(y, -1, 1)
            return x,y
    def wrap(self, data, lo, hi):
        data = n.where(data >= hi, lo + (data - hi), data)
        data = n.where(data < lo, hi + (data - lo), data)
        return data
    def drawmapboundary(self): pass
    def drawmeridians(self, locs, **kwargs):
        x,y = self(locs, n.zeros_like(locs))
        x = self.wrap(x, -2, 2)
        p.xticks(x, visible=False)
        p.grid(True)
    def drawparallels(self, lats, **kwargs):
        x,y = self(n.zeros_like(lats),lats)
        y = self.wrap(y, -1, 1)
        p.yticks(y, [str(L) for L in lats])
        p.grid(True)
    def makegrid(self, dim1, dim2, returnxy=True):
        y,x = n.indices((dim2,dim1))
        x = 4 * x.astype(n.float)/dim1 - 2
        y = 1 - 2 * y.astype(n.float)/dim2
        lon,lat = self(x,y, inverse=True)
        if returnxy: return lon,lat, x,y
        else: return lon,lat
    def imshow(self, *args, **kwargs):
        kwargs['extent'] = (-2,2,-1,1)
        return p.imshow(*args, **kwargs)

# Try to import basemap module, but on failure use the above class
try: from mpl_toolkits.basemap import Basemap
except(ImportError): 
    try: from matplotlib.toolkits.basemap import Basemap
    except(ImportError): pass

o = optparse.OptionParser()
o.set_usage('plot_map.py [options] mapfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cmap=True, max=True, drng=True)
o.add_option('-p', '--projection', dest='projection', default='moll',
    help='Map projection to use: moll (default), mill, cyl, robin, sinu.')
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plotting mode, can be log (default), lin.')
o.add_option('-c', '--cen', dest='cen', type='float', 
    help="Center longitude/right ascension (in degrees) of map.  Default is 0 for galactic coordinate output, 180 for equatorial.")
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
o.add_option('--osys', dest='osys', default='eq',
    help='Output coordinate system (plotted).  Can be eq (equatorial, default), ga (galactic), or ec (ecliptic)')
o.add_option('--iepoch', dest='iepoch', type='float', default=ephem.J2000,
    help='Epoch of input coordinates (in map).  Default J2000.')
o.add_option('--oepoch', dest='oepoch', type='float', default=ephem.J2000,
    help='Epoch of output coordinates (plotted).  Default J2000.')
o.add_option('--nobar', dest='nobar', action='store_true',
    help="Do not show colorbar.")
o.add_option('--res', dest='res', type='float', default=0.25,
    help="Resolution of plot (in degrees).  Default 0.25.")
o.add_option('--nside', dest='nside', type='int',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
o.add_option('--mask', dest='mask', type='float',
    help="Optional dB of weight below which data will be masked. Recommended=30")
o.add_option('--cat_gen',dest='cat',action='store_true',
    help="Optionally print out catalog of sources.")
opts,args = o.parse_args(sys.argv[1:])

cmap = p.get_cmap(opts.cmap)
if opts.cen is None:
    if opts.osys == 'eq': opts.cen = 180
    else: opts.cen = 0
map = Basemap(projection=opts.projection,lat_0=0,lon_0=opts.cen, rsphere=1.)
lons,lats,x,y = map.makegrid(360/opts.res,180/opts.res, returnxy=True)
# Mask off parts of the image to be plotted that are outside of the map
lt = lats[:,0]
ln1 = n.ones_like(lt) * (lons[lons.shape[0]/2,0])
ln2 = n.ones_like(lt) * (lons[lons.shape[0]/2,-1])
x1,y1 = map(ln1,lt); x2,y2 = map(ln2,lt)
x = n.ma.array(x)
for c,(i,j) in enumerate(zip(x1,x2)): x[c] = n.ma.masked_outside(x[c], i, j)
mask = x.mask
if opts.osys == 'eq': lons = 360 - lons
lats *= a.img.deg2rad; lons *= a.img.deg2rad
print 'Reading %s' % args[0]
h = a.map.Map(fromfits=args[0])
print 'SCHEME:', h.scheme()
print 'NSIDE:', h.nside()
if not opts.nside is None:
    nh = a.healpix.HealpixMap(nside=opts.nside)
    nh.from_hpm(h)
    h = nh
h.set_interpol(True)

crd = a.coord.radec2eq(n.array([lons.flatten(), lats.flatten()]))
m = a.coord.convert_m(opts.osys, opts.isys, 
    iepoch=opts.oepoch, oepoch=opts.iepoch)
x,y,z = n.dot(m, crd)
try: data, indices = h[x,y,z]
except(ValueError): data = h[x,y,z]
if not opts.mask is None:
    try:
        wgts = h.wgt[x,y,z]
        threshold = 10**(-opts.mask/10.)*n.max(wgts)
        msk = n.where(wgts > threshold, 1, 0)
        data *= msk
        print "Masking %2.0f%% of sky"% ((1 - msk.sum() / float(len(msk)))*100)
    except(AttributeError):
        print "Weights not included in file. No mask will be applied."
data.shape = lats.shape


# Generate source locations
if not opts.src is None:
    srclist,cutoff = a.scripting.parse_srcs(opts.src)
    cat = a.src.get_catalog(srcs=srclist, cutoff=cutoff)
    o = ephem.Observer()
    if opts.juldate is None:
        o.date = ephem.J2000
        o.epoch = o.date
        try: del(cat['Sun'])
        except(KeyError): pass
    else:
        o.date = a.phs.juldate2ephem(opts.juldate)
        o.epoch = o.date
    for s in cat.values():
        try: a.phs.RadioFixedBody.compute(s, o)
        except(TypeError): a.phs.RadioSpecial.compute(s, o)
    #cat.compute(o)
    # lat/lon coordinates of sources
    scrds = [ephem.Equatorial(s.ra,s.dec,epoch=o.epoch) for s in cat.values()]
    afreqs = n.array([.150])
    cat.update_jys(afreqs)
    sflxs = cat.get_jys().squeeze()
    snams = cat.keys()
    if opts.osys == 'ga':
        scrds = [ephem.Galactic(s, epoch=opts.oepoch) for s in scrds]
    elif opts.osys == 'ec':
        scrds = [ephem.Ecliptic(s, epoch=opts.oepoch) for s in scrds]
    slats = n.array([float(s.get()[1]) for s in scrds]) * a.img.rad2deg
    slons = n.array([float(s.get()[0]) for s in scrds]) * a.img.rad2deg
    if opts.osys == 'eq': slons = 360 - slons
    slons = n.where(slons < -180, slons + 360, slons)
    slons = n.where(slons >= 180, slons - 360, slons)

# Generate map grid/outline
map.drawmapboundary()
map.drawmeridians(n.arange(-180, 180, 30))
map.drawparallels(n.arange(-90,90,30)[1:], labels=[0,1,0,0], labelstyle='+/-')
# Set up data to plot
if opts.mode.startswith('log'): data = n.log10(n.abs(data))
if opts.max is None: max = data.max()
else: max = opts.max
if opts.drng is None:
    min = data.min()
    if min < (max - 10): min = max-10
else: min = max - opts.drng
data = data.clip(min, max)
data = n.ma.array(data, mask=mask)
map.imshow(data, vmax=max, vmin=min, cmap=cmap)


def mk_arr(val, dtype=n.double):
    if type(val) is n.ndarray: return val.astype(dtype)
    return n.array(val, dtype=dtype).flatten()
def get_val(x,y):
    lon,lat = map(x, y, inverse=True)
    if opts.osys == 'eq': lon = (360 - lon) % 360
    lon *= a.img.deg2rad; lat *= a.img.deg2rad
    ra,dec = ephem.hours(lon), ephem.degrees(lat)
    x,y,z = a.coord.radec2eq((ra,dec))
    #flx = h[(x,y,z)]
    crd = [mk_arr(c, dtype=n.double) for c in (x,y,z)]
    px,wgts = h.crd2px(*crd, **{'interpolate':1})
    flx = n.sum(h[px],axis=-1)
    return flx

# Plot src labels and markers on top of map image
if not opts.src is None:
    sx, sy = map(slons,slats)
    if not opts.cat is None: print "name, S_cat [Jy], S_avg [Jy]"
    for name, xpt, ypt, flx in zip(snams, sx, sy, sflxs):
        if xpt >= 1e30 or ypt >= 1e30: continue
        if opts.src_mark != '':
            map.plot(sx, sy, opts.src_color+opts.src_mark,markerfacecolor=None)
        if flx < 10: flx = 10
        p.text(xpt+.001, ypt+.001, name, size=5+2*int(n.round(n.log10(flx))),
            color=opts.src_color)
        if not opts.cat is None: 
            print "%s,%f,%f"%(name,flx,get_val(xpt,ypt))
if not opts.nobar: p.colorbar(shrink=.5, format='%.2f')
else: p.subplots_adjust(.05,.05,.95,.95)


if opts.outfile != '':
    print 'Saving to', opts.outfile
    p.savefig(opts.outfile)
else:
    # Add right-click functionality for finding locations/strengths in map.
    cnt = 1
    def click(event):
        global cnt
        if event.button == 3: 
            lon,lat = map(event.xdata, event.ydata, inverse=True)
            print event.xdata,type(event.xdata)
            if opts.osys == 'eq': lon = (360 - lon) % 360
            lon *= a.img.deg2rad; lat *= a.img.deg2rad
            ra,dec = ephem.hours(lon), ephem.degrees(lat)
            x,y,z = a.coord.radec2eq((ra,dec))
            flx = h[(x,y,z)]
            print '#%d (RA,DEC): (%s, %s), Jy: %f' % (cnt, ra, dec, flx)
            cnt += 1
        elif event.button==2:
            lon,lat = map(event.xdata, event.ydata, inverse=True)
            if opts.osys == 'eq': lon = (360 - lon) % 360
            lon *= a.img.deg2rad; lat *= a.img.deg2rad
            ra,dec = ephem.hours(lon), ephem.degrees(lat)
            x,y,z = a.coord.radec2eq((ra,dec))
            #flx = h[(x,y,z)]
            crd = [mk_arr(c, dtype=n.double) for c in (x,y,z)]
            px,wgts = h.crd2px(*crd, **{'interpolate':1})
            flx = n.sum(h[px],axis=-1)
            print '#%d (RA,DEC): (%s, %s), Jy: %f (4px sum)' % (cnt, ra, dec, flx)
            cnt += 1
        else: return
            


    #register this function with the event handler
    p.connect('button_press_event', click)
    p.show()
