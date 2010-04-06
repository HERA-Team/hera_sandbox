#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse, re
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
a.scripting.add_standard_options(o, max=True, drng=True)
opts, args = o.parse_args(sys.argv[1:])

RES = .01
SZ = (int(n.pi/2/RES), int(2*n.pi/RES)+2)
alt,az = n.indices(SZ)
alt = n.pi/2 - alt.astype(n.float) * RES
az = az.astype(n.float) * RES
tx,ty,tz = a.coord.azalt2top((az.flatten(), alt.flatten()))
map = Basemap(projection='ortho', lat_0=90, lon_0=180)
cx, cy = map(az * a.img.rad2deg, alt * a.img.rad2deg)

fname_regx = re.compile(r'.+_(\d+).fits')

m2 = int(n.sqrt(len(args)))
m1 = int(n.ceil(float(len(args)) / m2))
for cnt, filename in enumerate(args):
    print 'Reading', filename
    m = fname_regx.match(filename)
    bl = int(m.groups()[0])
    i,j = a.miriad.bl2ij(bl)
    h = a.map.Map(fromfits=filename)
    h.set_interpol(True)
    data = h[tx,ty,tz]
    data = 10*n.log10(n.abs(data))
    if not opts.max is None: max = opts.max
    else: max = data.max()
    if not opts.drng is None: min = max - opts.drng
    else: min = data.min()
    print max, min
    step = (max - min) / 20
    levels = n.arange(min-step, max+step, step)
    data.shape = az.shape
    p.subplot(m2, m1, cnt+1)
    map.drawmapboundary()
    map.drawmeridians(n.arange(0,360,30))
    map.drawparallels(n.arange(0,90,10))
    map.contourf(cx, cy, data.clip(min,max), levels)
    p.colorbar(format='%2.1f dB', shrink=.6)
    p.title(str((i,j)))

p.show()

