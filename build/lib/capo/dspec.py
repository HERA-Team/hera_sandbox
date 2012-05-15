import aipy as a, numpy as n
import capo as C

def sky_dly_thresh(bl_len, sdf, nchan, max_bl_frac=1.):
    '''For bl_len in ns, return the (upper,lower) delay bins that geometrically
    correspond to the sky.'''
    bin_dly = 1. / (sdf * nchan)
    max_bl = bl_len * max_bl_frac
    uthresh, lthresh = max_bl/bin_dly + 1.5, -max_bl/bin_dly - 0.5
    uthresh, lthresh = int(n.ceil(uthresh)), int(n.floor(lthresh))
    return (uthresh,lthresh)
    
def all_sky_dly_thresh(aa, sdf, nchan, max_bl_frac=1.):
    '''Return a dictionary, indexed by baseline, of the (upper,lower) delay 
    bins that geometrically correspond to the sky.'''
    filters = {}
    for i in range(len(aa.ants)):
      for j in range(len(aa.ants)):
        if j < i: continue
        bl = aa.ij2bl(i,j)
        max_bl = aa.get_baseline(i,j)
        max_bl = n.sqrt(n.dot(max_bl, max_bl))
        filters[bl] = sky_dly_thresh(max_bl, sdf, nchan, max_bl_frac)
    return filters

def wideband_dspec(jy_spec, wgts, uthresh, lthresh, tol=1e-4, window='none'):
    window = a.dsp.gen_window(jy_spec.size, window=window)
    _d = n.fft.ifft(jy_spec * window)
    _w = n.fft.ifft(wgts * window)
    area = n.ones(_d.size, dtype=n.int); area[uthresh:lthresh] = 0
    _d_cl, info = a.deconv.clean(_d, _w, area=area, tol=tol, stop_if_div=False)
    d_mdl = n.fft.fft(_d_cl)
    d_res = jy_spec - d_mdl * wgts
    if False:
        import pylab as p
        p.semilogy(n.abs(_d), 'g')
        p.semilogy(n.abs(_d_cl), 'k.')
        p.semilogy(n.abs(n.fft.ifft(d_res * window)), 'r')
        p.show()
    return d_mdl, d_res

