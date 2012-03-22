#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys

if not 'autos.npz' in sys.argv:
    data, wgts = {}, {}
    for filename in sys.argv[1:]:
        print filename
        t,d,f = C.arp.get_dict_of_uv_data([filename], 'auto', 'xx')
        for bl in d:
            data[bl] = data.get(bl,[]) + [d[bl].sum(axis=0)]
            wgts[bl] = wgts.get(bl,[]) + [n.logical_not(f[bl]).sum(axis=0)]
    data = n.array(data.values(), dtype=n.float32)
    wgts = n.array(wgts.values(), dtype=n.int16)

    n.savez('autos.npz', data=data, wgts=wgts)

d = n.load('autos.npz')
data, wgts = d['data'], d['wgts']
d = n.where(wgts > wgts.max()/2, data/wgts, 0)
w = n.where(wgts > wgts.max()/2, 1., 0)
sdf = 4.8828125e-5
freqs = .1 + n.arange(data.shape[-1]) * sdf

if True: # Global Signal mode
    data *= 30. # approximate gain cal in K
    if True:
        d,w = 0,0
        for i in range(data.shape[0]):
            #if i in [40,55]: continue
            if not i in [1]: continue
            d += data[i].sum(axis=0)
            w += wgts[i].sum(axis=0)
        _d, _w = n.fft.ifft(d), n.fft.ifft(w)
        _c,info = a.deconv.clean(_d,_w, tol=0, verbose=True, maxiter=30000)
        p.semilogy(n.abs(_c))
        #p.plot(freqs, n.fft.fft(_c))
        #p.plot(freqs, d / n.where(w==0,1,w) - n.fft.fft(_c))
    else:
        for i in range(data.shape[0]):
            d = data[i].sum(axis=0)
            w = wgts[i].sum(axis=0)
            p.plot(freqs, d / n.where(w == 0, 1, w))
    p.show()
    sys.exit(0)

_d = []
_d_sum, _d_wgt = 0, 0
for i in range(d.shape[0]):
#for i in range(4):
    if i in [40,55]: continue
    d[i] *= 30. # approximate gain cal in K
    print i, d[i].shape, n.average(d[i])
    x = n.arange(d[i].shape[0])
    y = n.arange(d[i].shape[1])
    if True:
        di_smooth = []
        for t in range(d[i].shape[0]):
            d_ = d[i,t,300:-250]
            w_ = w[i,t,300:-250]
            y_ = y[300:-250]
            if n.all(w_ == 0): continue
            poly = n.polyfit(y_.compress(w_), n.log10(d_.compress(w_)), deg=6)
            d[i,t,:] -= 10**n.polyval(poly, y) * w[i,t,:]
            di_smooth.append(10**n.polyval(poly,y))
        di_smooth = n.array(di_smooth)[...,300:-250]
    if True:
        for c in range(d[i].shape[1]):
            d_ = d[i,:,c]
            w_ = w[i,:,c]
            if n.all(w_ == 0): continue
            if False: # remove time average
                d[i,:,c] -= d_.sum() / w_.sum() * w_
            elif False: # remove poly in time
                poly = n.polyfit(x.compress(w_), d_.compress(w_), deg=6)
                d[i,:,c] -= n.polyval(poly, x) * w[i,:,c]
            else:
                d[i,1:-1,c] -= .5 * (d[i,:-2,c] + d[i,2:,c])
                w[i,0,c] = 0
                w[i,-1,c] = 0
                w[i,1:-1,c] = w[i,:-2,c] * w[i,2:,c]
                d[i,:,c] *= w[i,:,c]
    di = d[i,:,300:-250]
    wi = w[i,:,300:-250]
    # Flag some more
    std = n.sqrt(n.sum(di**2)/n.sum(wi))
    flagi = n.where(n.abs(di) > 3*std, 0, 1)
    di *= flagi
    wi *= flagi
    # Think about whitening spectrum here
    #di = n.where(wi == 0, n.random.normal(scale=di_smooth/n.sqrt(10.7*4.88e4)), di)
    if False:
        C.arp.waterfall(n.where(wi, di, -20), mode='real', mx=5, drng=10)
        p.colorbar()
        p.show()
        sys.exit(0)
    if False: # Use whole-band FFT
        window = a.dsp.gen_window(di.shape[-1], window='kaiser3')
        _di = n.fft.ifft(di * window)
        _wi = n.fft.ifft(wi * window)
        gi = n.sqrt(n.average((wi*window)**2, axis=1))
        for t in range(_di.shape[0]):
            _d_c,info = a.deconv.clean(_di[t], _wi[t], tol=1e-4)
            _di[t] = _d_c + info['res'] / gi[t]
        B = freqs[-250] - freqs[300]
        fq = n.average(freqs[300:-250])
    else: # Use PFB analysis
        SUBBAND = 300
        ch_cen = 1000
        taps = 3
        ch0,ch1 = ch_cen - taps*SUBBAND/2, ch_cen + taps*SUBBAND/2
        _di = C.pfb.pfb(di[:,ch0:ch1],taps=taps,window='kaiser3',fft=n.fft.ifft)
        _wi = C.pfb.pfb(wi[:,ch0:ch1],taps=taps,window='kaiser3',fft=n.fft.ifft)
        if True: # Deconv by sampling pattern
            gi = n.sqrt(n.average((wi[:,ch_cen-SUBBAND/2:ch_cen+SUBBAND/2])**2, axis=1))
            for t in range(_di.shape[0]):
                _d_c,info = a.deconv.clean(_di[t], _wi[t], tol=1e-4)
                _di[t] = _d_c + info['res'] / gi[t]
        B = (freqs[300+ch1] - freqs[300+ch0]) / taps
        fq = n.average(freqs[300+ch0:300+ch1])
    for _dj in _d:
        _d_sum += _di * n.conj(_dj)
        _d_wgt += 1
    _d.append(_di)
    if False:
        C.arp.waterfall(a.img.recenter(_di, (0,_di.shape[-1]/2)), mx=-1, drng=2)
        p.colorbar()
        p.show()
        sys.exit(0)

print B, fq

_d = _d_sum / _d_wgt

if True:
    C.arp.waterfall(a.img.recenter(_d, (0,_d.shape[-1]/2)), mx=-3, drng=3)
    p.colorbar()
    p.show()
    sys.exit(0)

#good_ints = range(22,31) + range(32,38)
good_ints = range(25,37)

_d_sum = 0
for i in good_ints: _d_sum += _d[i]
_d = _d_sum / len(good_ints)

# some cosmology
z = C.pspec.f2z(fq)
eta = n.fft.fftfreq(_d.size, sdf)
eta = eta[:_d.size/2]
k = C.pspec.dk_deta(z) * eta
_d = _d[:_d.size/2]
k = n.sqrt(k**2 + .006**2)
scalar = C.pspec.k3pk_from_Trms(1, k=k, fq=fq, B=B)
_d *= scalar
p.loglog(k, n.abs(_d), 'b.')
k,_d = C.pspec.rebin_log(k, _d, nbins=10)
p.loglog(k, n.abs(_d), 'k-')
p.show()

