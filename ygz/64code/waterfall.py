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
def waterfall(d, ax=None, mode='log', mx=None, drng=None, recenter=False, **kwargs):
    if n.ma.isMaskedArray(d): d = d.filled(0)
    if recenter: d = a.img.recenter(d, n.array(d.shape)/2)
    d = data_mode(d, mode=mode)
    if mx is None: mx = d.max()
    if drng is None: drng = mx - d.min()
    mn = mx - drng
    if ax:
    	return ax.imshow(d, vmax=mx, vmin=mn, aspect='auto', interpolation='nearest', **kwargs)
    else:
    	return plt.imshow(d, vmax=mx, vmin=mn, aspect='auto', interpolation='nearest', **kwargs)