#!/usr/bin/env python

import aipy as a, numpy as n, pylab as p, sys, os, ephem, optparse
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_beam.py [options] mapfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cmap=True, max=True, drng=True,cal=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('--res', dest='res', type='float', default=0.25,
    help="Resolution of plot (in degrees).  Default 0.25.")
o.add_option('-i', '--interp',  action='store_true', 
    help="Interpolate")

def data_mode(data, mode='abs'):
    if mode.startswith('phs'): data = n.angle(data)
    elif mode.startswith('lin'):
        data = n.absolute(data)
        data = n.masked_less_equal(data, 0)
    elif mode.startswith('real'): data = data.real
    elif mode.startswith('imag'): data = data.imag
    elif mode.startswith('log'):
        data = n.absolute(data)
        data = n.log10(data)
    else: raise ValueError('Unrecognized plot mode.')
    return data

opts,args = o.parse_args(sys.argv[1:])

cmap = p.get_cmap(opts.cmap)

map = Basemap(projection='ortho',lat_0=90,lon_0=180,rsphere=1.)
h = a.map.Map(fromfits=args[0])
print 'SCHEME:', h.scheme()
print 'NSIDE:', h.nside()
if opts.interp: h.set_interpol(True)

lons,lats,x,y = map.makegrid(360/opts.res,180/opts.res, returnxy=True)
lons = 360 - lons
lats *= a.img.deg2rad; lons *= a.img.deg2rad
y,x,z = a.coord.radec2eq(n.array([lons.flatten(), lats.flatten()]))
ax,ay,az = a.coord.latlong2xyz(n.array([0,0]))
#x = x - ax
#y = y - ay
#z = z - az
try: data, indices = h[x,y,z]
except(ValueError): data = h[x,y,z]
data.shape = lats.shape
data /= h[0,0,1]
#data = data**2 # only if a voltage beam
data = data_mode(data, opts.mode)


map.drawmapboundary()
map.drawmeridians(n.arange(0, 360, 30))
map.drawparallels(n.arange(0, 90, 10))

if opts.max is None: max = data.max()
else: max = opts.max
if opts.drng is None:
    min = data.min()
#    if min < (max - 10): min = max-10
else: min = max - opts.drng
step = (max - min) / 10
levels = n.arange(min-step, max+step, step)
print min,max
#data = data.clip(min, max)
#data = n.ma.array(data, mask=mask)
#min=0
map.imshow(data, vmax=max, vmin=min, cmap=cmap)
#map.contourf(cx,cy,data,levels,linewidth=0,cmap=cmap)
p.colorbar()

p.show()
