import aipy as a, numpy as n, pylab as plt
import sys, scipy
from mpl_toolkits.basemap import Basemap

def data_mode(data, mode='abs'):
    if mode.startswith('phs'): data = n.angle(data)
    elif mode.startswith('lin'):
        data = n.absolute(data)
        #data = n.masked_less_equal(data, 0)
    elif mode.startswith('real'): data = data.real
    elif mode.startswith('imag'): data = data.imag
    elif mode.startswith('log'):
        data = n.absolute(data)
        data = n.log10(data)
    else: raise ValueError('Unrecognized plot mode.')
    return data

def waterfall(d, mode='log', mx=None, drng=None, recenter=False, **kwargs):
    if n.ma.isMaskedArray(d): d = d.filled(0)
    if recenter: d = a.img.recenter(d, n.array(d.shape)/2)
    d = data_mode(d, mode=mode)
    if mx is None: mx = d.max()
    if drng is None: drng = mx - d.min()
    mn = mx - drng
    return plt.imshow(d, vmax=mx, vmin=mn, aspect='auto', interpolation='nearest', **kwargs)

def plot_hmap_ortho(h, cmap='jet', mode='log', mx=None, drng=None, 
        res=0.25, verbose=False, normalize=False):
    m = Basemap(projection='ortho',lat_0=90,lon_0=180,rsphere=1.)
    if verbose:
        print 'SCHEME:', h.scheme()
        print 'NSIDE:', h.nside()
    lons,lats,x,y = m.makegrid(360/res,180/res, returnxy=True)
    lons = 360 - lons
    lats *= a.img.deg2rad; lons *= a.img.deg2rad
    y,x,z = a.coord.radec2eq(n.array([lons.flatten(), lats.flatten()]))
    ax,ay,az = a.coord.latlong2xyz(n.array([0,0]))
    data = h[x,y,z]
    data.shape = lats.shape
    if normalize: data /= h[0,0,1]
    data = data_mode(data, mode)
    m.drawmapboundary()
    m.drawmeridians(n.arange(0, 360, 30))
    m.drawparallels(n.arange(0, 90, 10))
    if mx is None: mx = data.max()
    if drng is None:
        mn = data.min()
    #    if min < (max - 10): min = max-10
    else: mn = mx - drng
    return m.imshow(data, vmax=mx, vmin=mn, cmap=cmap)
    #step = (mx - mn) / 10
    #levels = n.arange(mn-step, mx+step, step)
    #map.contourf(cx,cy,data,levels,linewidth=0,cmap=cmap)
    
