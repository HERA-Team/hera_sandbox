import aipy as a, numpy as np, pylab as plt
import sys, scipy
from mpl_toolkits.basemap import Basemap

def data_mode(data, mode='abs'):
    if mode.startswith('phs'): data = np.angle(data)
    elif mode.startswith('lin'):
        data = np.absolute(data)
        #data = np.masked_less_equal(data, 0)
    elif mode.startswith('real'): data = data.real
    elif mode.startswith('imag'): data = data.imag
    elif mode.startswith('log'):
        data = np.absolute(data)
        data = np.log10(data)
    else: raise ValueError('Unrecognized plot mode.')
    return data

def waterfall(d, mode='log', mx=None, drng=None, recenter=False, **kwargs):
    if np.ma.isMaskedArray(d): d = d.filled(0)
    if recenter: d = a.img.recenter(d, np.array(d.shape)/2)
    d = data_mode(d, mode=mode)
    if mx is None: mx = d.max()
    if drng is None: drng = mx - d.min()
    mn = mx - drng
    return plt.imshow(d, vmax=mx, vmin=mn, aspect='auto', interpolation='nearest', **kwargs)

def plot_hmap_ortho(h, cmap='jet', mode='log', mx=None, drng=None, 
        res=0.25, verbose=False):
    m = Basemap(projection='ortho',lat_0=90,lon_0=180,rsphere=1.)
    if verbose:
        print 'SCHEME:', h.scheme()
        print 'NSIDE:', h.nside()
    lons,lats,x,y = m.makegrid(360/res,180/res, returnxy=True)
    lons = 360 - lons
    lats *= a.img.deg2rad; lons *= a.img.deg2rad
    y,x,z = a.coord.radec2eq(np.array([lons.flatten(), lats.flatten()]))
    ax,ay,az = a.coord.latlong2xyz(np.array([0,0]))
    data = h[x,y,z]
    data.shape = lats.shape
    data /= h[0,0,1]
    data = data_mode(data, mode)
    m.drawmapboundary()
    m.drawmeridians(np.arange(0, 360, 30))
    m.drawparallels(np.arange(0, 90, 10))
    if mx is None: mx = data.max()
    if drng is None:
        mn = data.min()
    #    if min < (max - 10): min = max-10
    else: mn = mx - drng
    return m.imshow(data, vmax=mx, vmin=mn, cmap=cmap)
    #step = (mx - mn) / 10
    #levels = np.arange(mn-step, mx+step, step)
    #map.contourf(cx,cy,data,levels,linewidth=0,cmap=cmap)
    
def omni_view(reds,vis,pol,int=10,chan=500,norm=False,cursor=True,save=None,colors=None,symbols=None):
    if not colors:
        colors = ["#006BA4", "#FF7F0E", "#2CA02C", "#D61D28", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"]
    if not symbols: 
        symbols = ["o", "v", "^", "<", ">", "*"]
    points = []
    sym = []
    col = []
    bl = []
    ngps = len(reds)
    if save:
        p.clf()
        p.cla()
    for i,gp in enumerate(reds):
        c = colors[i%len(colors)]
        s = symbols[i/len(colors)]
        for r in gp:
            try:
                points.append(vis[r][pol][int,chan])
                bl.append(r)
            except(KeyError):
                points.append(np.conj(vis[r[::-1]][pol][int,chan]))
                bl.append(r[::-1])
            sym.append(s)
            col.append(c)
    points = np.array(points)
    max_x=0
    max_y=0
    for i,pt in enumerate(points):
        if norm:
            p.scatter(pt.real/np.abs(pt), pt.imag/np.abs(pt), c=col[i], marker=sym[i], s=50, label='{}'.format(bl[i]))
        else:        
            p.scatter(pt.real, pt.imag, c=col[i], marker=sym[i], s=50, label='{}'.format(bl[i]))
            if np.abs(pt.real) > max_x: max_x = np.abs(pt.real)
            if np.abs(pt.imag) > max_y: max_y = np.abs(pt.imag)
            
    if norm:         
        p.xlim(-1,1)
        p.ylim(-1,1)
    else: 
        p.xlim(-max_x-.1*max_x,max_x+.1*max_x)
        p.ylim(-max_y-.1*max_y,max_y+.1*max_y)
    p.ylabel('imag(V)')
    p.xlabel('real(V)')
    if cursor:
        from mpldatacursor import datacursor
        datacursor(formatter='{label}'.format)
    if save:
        p.savefig(save)
    return None

def omni_view_gif(filenames, name='omni_movie.gif'):
    '''Give it a full path for images and all the images you want to use in order to use. list format.'''
    import imageio
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))
        imageio.mimsave(name, images)
