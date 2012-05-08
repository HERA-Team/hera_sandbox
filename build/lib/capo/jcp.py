import os, numpy as n, pylab as p, aipy as a

def read_srcnpz(npzfiles, verbose=False, with_xyz=True):
    times, fluxes, wgts, x,y,z = {}, {}, {}, {},{},{}
    for npz in npzfiles:
        if verbose: print 'Reading', npz
        f = open(npz)
        src = os.path.basename(npz).split('__')[0]
        try: npz = n.load(f)
        except:
            print 'Load file failed'
            continue
        try: wgt = npz['wgts']
        except(KeyError):
            try: wgt = n.logical_not(npz['mask'])
            except(KeyError): wgt = n.ones_like(npz['times'])
        times[src] = times.get(src,[]) + [npz['times']]
        fluxes[src] = fluxes.get(src,[]) + [npz['spec']]
        wgts[src] = wgts.get(src,[]) + [wgt]
        if with_xyz:
            x[src] = x.get(src,[]) + [npz['x']]
            y[src] = y.get(src,[]) + [npz['y']]
            z[src] = z.get(src,[]) + [npz['z']]
        f.close()
    for src in times:
        times[src] = n.concatenate(times[src])
        fluxes[src] = n.concatenate(fluxes[src])
        wgts[src] = n.concatenate(wgts[src])
        if with_xyz:
            x[src] = n.concatenate(x[src])
            y[src] = n.concatenate(y[src])
            z[src] = n.concatenate(z[src])
    if with_xyz: return times, fluxes, wgts, x, y, z
    else: return times, fluxes, wgts

def plot_hpx(hpx, mode='log', pol='x', mx=None, drng=2, interpolation=False, 
        colorbar=False, nogrid=False):
    im = a.img.Img(size=400, res=.4)
    x,y,z = im.get_top(center=n.array(im.shape)/2)
    assert(pol in 'xy')
    hpx.set_interpol(interpolation)
    if pol == 'x': d = hpx[x,y,z]
    else: d = hpx[-y,x,z]
    d.shape = x.shape
    d = n.where(x.mask, 0, d)
    if mode.startswith('phs'): d = n.angle(d)
    elif mode.startswith('lin'): d = n.absolute(d)
    elif mode.startswith('real'): d = d.real
    elif mode.startswith('imag'): d = d.imag
    elif mode.startswith('log'):
        d = n.absolute(d)
        d = n.ma.masked_less_equal(d, 0)
        d = n.ma.log10(d)
    else: raise ValueError('Unrecognized plot mode.')
    if mx is None: mx = d.max()
    if drng is None: drng = mx - d.min()
    mn = mx - drng
    if not nogrid:
        from mpl_toolkits.basemap import Basemap
        map = Basemap(projection='ortho', lon_0=180, lat_0=90,
            rsphere=1, llcrnrx=-1.25, llcrnry=-1.25, urcrnrx=1.25,urcrnry=1.25)
        map.drawmeridians(n.arange(-180,180,30))
        map.drawparallels(n.arange(0,90,15))
        map.drawmapboundary()
        map.imshow(d, vmin=mn, vmax=mx, interpolation='nearest')
    else: p.imshow(d, vmin=mn, vmax=mx, origin='lower', interpolation='nearest')

    #data.shape = az.shape
    #m.drawmapboundary()
    #m.drawmeridians(n.arange(0,360,30))
    #m.drawparallels(n.arange(0,90,10))
    #step = (max - min) / 35.
    #levels = n.arange(min-step, max+step, step)
    #m.contourf(cx, cy, data, levels)
    if colorbar: p.colorbar(shrink=.5)
    #return m

