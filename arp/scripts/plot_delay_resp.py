#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p
from mpl_toolkits.basemap import Basemap

FQ = .164
BINWIDTH = 1. / 7200. # Hz
bin_edges = n.arange(-200,200,1)
NBINS = len(bin_edges)

def mk_dly(bl, ex, ey, ez):
    #f = -(bl[0]*n.sin(lons) + bl[1]*n.cos(lons)) * n.cos(lats)
    return (bl[0]*ex + bl[1]*ey + bl[2]*ez)

aa = a.cal.get_aa('psa898_v003', n.array([FQ]))
if False:
    for ant in aa: ant.set_pointing(twist=n.pi/4)
if False: aa.lat = '0:00'; aa.update()
lat = aa.lat
print lat
fig = p.figure(figsize=(4,5))
fig.subplots_adjust(left=.05, top=.95, bottom=.05, right=0.95)
m = Basemap(projection='ortho', lat_0=lat, lon_0=180, rsphere=1.)
h = a.healpix.HealpixMap(nside=128)
DIM = 400
RES = .5
SIZE = DIM / RES
im = a.img.Img(DIM, res=RES)

def proj_hmap(dat, mode='e'):
    h.map = dat
    #if mode == 'e': x,y,z = im.get_eq(dec=aa.lat, center=(SIZE/2,SIZE/2))
    if mode == 'e': x,y,z = im.get_eq(dec=lat, center=(SIZE/2,SIZE/2))
    elif mode == 't': x,y,z = im.get_top(center=(SIZE/2,SIZE/2))
    else: raise ValueError(mode)
    v = n.logical_not(x.mask)
    x,y,z = x.flatten(), y.flatten(), z.flatten()
    d = h[x,y,z]
    d.shape = (SIZE,SIZE)
    d = n.where(v, d, n.NaN)
    return d

xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
tx,ty,tz = n.dot(aa._eq2zen, xyz)
_bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]
_bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]
bm = 0.5 * (_bmx**2 + _bmy**2)
ones = n.where(tz > 0, 1, 0)
bm = n.where(tz > 0, bm, 0)
bl = aa.get_baseline(0,16,'r')
dly = mk_dly(bl, *xyz)
#hist, bin_edges = n.histogram(dly, bins=bin_edges, weights=bm**2)
hist, bin_edges = n.histogram(dly, bins=bin_edges, weights=ones**2)
bins = 0.5 * (bin_edges[1:] + bin_edges[:-1])

p.plot(bins, hist); p.show()

m.imshow(proj_hmap(dly,'e'), origin='lower')

m.drawmapboundary(linewidth=2)
m.drawmeridians(n.arange(0, 360, 30), linewidth=2)
m.drawparallels(n.arange(-90,90,30), linewidth=2)
m.colorbar()

p.show()
