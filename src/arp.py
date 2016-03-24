import aipy as a, numpy as n, pylab as plt
import sys, scipy
from plot import data_mode, waterfall, plot_hmap_ortho

def get_dict_of_uv_data(filenames, antstr, polstr, decimate=1, decphs=0, verbose=False, recast_as_array=True):
    info = {'lsts':[], 'times':[]}
    ts = {}
    dat, flg = {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        if decimate > 1: uv.select('decimate', decimate, decphs)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            if not ts.has_key(t):
                info['times'].append(t)
                info['lsts'].append(uv['lst'])
                ts[t] = None
            #bl = a.miriad.ij2bl(i,j)
            bl = (i,j)
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            pol = a.miriad.pol2str[uv['pol']]
            if not dat[bl].has_key(pol):
                dat[bl][pol],flg[bl][pol] = [],[]
            dat[bl][pol].append(d)
            flg[bl][pol].append(f)
    if recast_as_array:
        # This option helps reduce memory footprint, but it shouldn't
        # be necessary: the replace below should free RAM as quickly
        # as it is allocated.  Unfortunately, it doesn't seem to...
        for bl in dat.keys():
          for pol in dat[bl].keys():
            dat[bl][pol] = n.array(dat[bl][pol])
            flg[bl][pol] = n.array(flg[bl][pol])
        info['lsts'] = n.array(info['lsts'])
        info['times'] = n.array(info['times'])
    return info, dat, flg

def clean_transform(d, w=None, f=None, clean=1e-3, window='blackman-harris'):
    #d = d.swapaxes(0, axis)
    #f = n.logical_not(f.swapaxes(0, axis))
    if w is None and not f is None: w = n.logical_not(f)
    elif w is None: w = n.ones(d.shape, dtype=n.float)
    window = a.dsp.gen_window(d.shape[-1], window=window)
    _d = n.fft.ifft(d*window, axis=-1)
    _w = n.fft.ifft(w*window, axis=-1)
    if _d.ndim == 2:
        for i in range(_d.shape[0]): # XXX would be nice to make this work on any shape of d
            g = n.sqrt(n.average(w[i]**2))
            if g == 0: continue
            _d[i],info = a.deconv.clean(_d[i], _w[i], tol=clean)
            _d[i] += info['res'] / g
    else: 
        g = n.sqrt(n.average(w**2))
        if g != 0:
            _d,info = a.deconv.clean(_d, _w, tol=clean)
            _d += info['res'] / g
    #_d = _d.swapaxes(0, axis)
    return _d

def gen_ddr_filter(shape, dw, drw, ratio=.25, invert=False):
    filter = n.ones(shape)
    x1,x2 = drw, -drw
    if x2 == 0: x2 = shape[0]
    y1,y2 = dw, -dw
    if y2 == 0: y2 = shape[1]
    filter[x1+1:x2,0] = 0
    filter[0,y1+1:y2] = 0
    filter[1:,1:] = 0
    x,y = n.indices(shape).astype(n.float)
    x -= shape[0]/2
    y -= shape[1]/2
    r2 = (x/(ratio*drw+.5))**2 + (y/(ratio*dw+.5))**2
    r2 = a.img.recenter(r2, (shape[0]/2, shape[1]/2))
    filter += n.where(r2 <= 1, 1, 0)
    filter = filter.clip(0,1)
    if invert: return n.logical_not(filter)
    else: return filter

def rms(d,wgt=None):
    if wgt == None: return n.sqrt(n.average(n.abs(d)**2))
    else: return n.sqrt(n.sum(n.abs(d)**2) / n.sum(n.abs(wgt)**2))

def redundant_bl_cal(d1, w1, d2, w2, fqs, use_offset=False, maxiter=10, window='blackman-harris',
        clean=1e-4, verbose=False, tau=0., off=0.):
    '''Return gain and phase difference between two redundant measurements
    d1,d2 with respective weights w1,w2.'''
    # Compute measured values
    dtau,doff,mx = 0,0,0
    d12 = d2 * n.conj(d1)
    # For 2D arrays, assume first axis is time and integrate over it
    if d12.ndim > 1: d12_sum,d12_wgt = n.sum(d12,axis=0), n.sum(w1*w2,axis=0)
    else: d12_sum,d12_wgt = d12, w1*w2
    if n.all(d12_wgt == 0): return n.zeros_like(d12_sum), 0.
    d11 = d1 * n.conj(d1)
    if d11.ndim > 1: d11_sum,d11_wgt = n.sum(d11,axis=0), n.sum(w1*w1,axis=0)
    else: d11_sum,d11_wgt = d11, w1*w1
    window = a.dsp.gen_window(d12_sum.size, window=window)
    dlys = n.fft.fftfreq(fqs.size, fqs[1]-fqs[0])
    # Begin at the beginning
    d12_sum *= n.exp(-2j*n.pi*(fqs*tau+off))
    for j in range(maxiter):
        d12_sum *= n.exp(-2j*n.pi*(fqs*dtau+doff))
        tau += dtau; off += doff
        _phs = n.fft.fft(window*d12_sum)
        _wgt = n.fft.fft(window*d12_wgt)
        _phs,info = a.deconv.clean(_phs, _wgt, tol=clean)
        #_phs += info['res'] / a.img.beam_gain(_wgt)
        _phs = n.abs(_phs)
        mx = n.argmax(_phs)
        if j > maxiter/2 and mx == 0: # Fine-tune calibration with linear fit
            valid = n.where(d12_wgt > d12_wgt.max()/2, 1, 0)
            valid *= n.where(n.abs(d12_sum) > 0, 1, 0) # Throw out zeros, which NaN in the log below
            fqs_val = fqs.compress(valid)
            dly = n.real(n.log(d12_sum.compress(valid))/(2j*n.pi)) # This doesn't weight data
            wgt = d12_wgt.compress(valid); wgt.shape = (wgt.size,1)
            B = n.zeros((fqs_val.size,1)); B[:,0] = dly
            if use_offset: # allow for an offset component
                A = n.zeros((fqs_val.size,2)); A[:,0] = fqs_val; A[:,1] = 1
                dtau,doff = n.linalg.lstsq(A*wgt,B*wgt)[0].flatten()
            else:
                #A = n.zeros((fqs_val.size,1)); A[:,0] = fqs_val
                #dtau = n.linalg.lstsq(A*wgt,B*wgt)[0].flatten()[0]
                dtau = n.sum(wgt.flatten()*dly/fqs_val) / n.sum(wgt.flatten())
        else: # Pull out an integral number of phase wraps
            if mx > _phs.size/2: mx -= _phs.size
            dtau,doff = mx / (fqs[-1] - fqs[0]), 0
            mxs = mx + n.array([-1,0,1])
            dtau = n.sum(_phs[mxs] * dlys[mxs]) / n.sum(_phs[mxs])
            #dtau = n.sum(_phs**2 * dlys) / n.sum(_phs**2)
            #dtau = n.sum(_phs * dlys) / n.sum(_phs)
        if verbose: print j, dtau, doff, (tau, off), mx
        #P.subplot(211); P.plot(n.fft.fftshift(dlys), n.fft.fftshift(_phs)); P.xlim(-200,200)
        #P.subplot(212); P.plot(fqs, n.angle(d12_sum))
        #P.show()
    #P.subplot(211); P.plot(n.fft.fftshift(dlys), n.fft.fftshift(_phs)); P.xlim(-200,200)
    #P.subplot(212); P.plot(fqs, n.angle(d12_sum))
    #P.show()
    off %= 1
    info = {'dtau':dtau, 'doff':doff, 'mx':mx} # Some information about last step, useful for detecting screwups
    g12 = d12_sum / d12_wgt.clip(1,n.Inf)
    g11 = d11_sum / d11_wgt.clip(1,n.Inf)
    gain = n.where(g11 != 0, g12/g11, 0)
    if use_offset: return gain, (tau,off), info
    else: return gain, tau, info

def selfcal_diff(m_bl, r_ant, r_wgt=1e6):
    ants = {}
    for bl in m_bl:
        i,j = a.miriad.bl2ij(bl)
        ants[i] = ants[j] = None
    ants = ants.keys(); ants.sort()
    def antind(ant): return ants.index(ant)
    NANT = len(ants)
    NMEAS = len(m_bl) + len(r_ant)
    P = n.zeros((NMEAS, NANT), dtype=n.double)
    M = n.zeros((NMEAS, 1), dtype=n.double)
    # Put in reference information (i.e. assumptions)
    for cnt1,(i,val) in enumerate(r_ant.items()):
        P[cnt1,antind(i)] = r_wgt
        M[cnt1,0] = r_wgt * val
    cnt1 += 1 # Set cnt1 to point to the next slot
    # Add in measurements
    for cnt2,(bl,val) in enumerate(m_bl.items()):
        i,j = a.miriad.bl2ij(bl)
        P[cnt1+cnt2,antind(j)] = 1; P[cnt1+cnt2,antind(i)] = -1
        M[cnt1+cnt2,0] = val
    # Now that information is in matrix form, solve it
    Pinv = n.linalg.pinv(P) # this succeeds where lstsq fails, for some reason
    C = n.dot(Pinv, M)
    return dict(zip(ants,C))

def sinuspike(d, fqs, f=None, clean=1e-3, maxiter=100, nsig=3, window='blackman-harris'):
    d = d * (fqs/.150)**2.5 # Flatten noise assuming synchrotron spectral index
    sig = n.sqrt(n.median(n.abs(d)**2))
    if sig == 0: return
    _mdl = n.zeros_like(d)
    if not f is None: w = n.logical_not(f).astype(n.float)
    else: w = n.ones(d.shape, dtype=n.float)
    for i in range(maxiter):
        di = d - w*n.fft.fft(_mdl)
        sig = n.sqrt(n.median(n.abs(di)**2))
        #_di = clean_transform(di, w=w, clean=clean, window=window)
        _di = n.fft.ifft(di)
        _sig = n.sqrt(n.median(n.abs(_di)**2))
        adi,_adi = n.abs(di), n.abs(_di)
        amx, _amx = n.argmax(adi), n.argmax(_adi)
        mx,_mx = adi[amx]/sig, _adi[_amx]/_sig
        if mx < nsig and _mx < nsig: break
        if mx > _mx: d[amx],w[amx] = 0, 0
        else: _mdl[_amx] += .3 * _di[_amx]
    f = n.logical_not(w).astype(n.int)
    mdl = n.fft.fft(_mdl) * (fqs/.150)**-2.5
    return mdl, f
