#! /usr/bin/env python
import aipy as a, numpy as n, pylab as P
import capo as C
import sys

uv = a.miriad.UV(sys.argv[-1])
sdf = uv['sdf']
freqs = uv['sfreq'] + n.arange(uv['nchan'])*sdf
del(uv)
t,dat,flg = C.arp.get_dict_of_uv_data(sys.argv[1:], 'auto', 'xx')
#nants = len(ants.split(','))
#nants = 64
nants = 9

ncols = n.floor(n.sqrt(nants))
nrows = n.ceil(n.sqrt(nants))

for i,bl in enumerate(dat):
    print bl
    d = n.ma.array(dat[bl],mask=flg[bl])
    d = n.ma.mean(d, axis=0)
    d_df = d.copy(); d_df[1:-1] -= (d_df[:-2] + d_df[2:])/2
    sig = n.sqrt(n.median(n.abs(d_df)**2))
    #print sig
    d = n.where(n.abs(d_df) > 10*sig, 0, d)
    d = d[300:-250]
    d_df = d_df[300:-250]
    cfreqs = freqs[300:-250]
    w = n.where(d == 0, 0, 1.)
    if True:
        poly = n.polyfit(cfreqs.compress(w), n.log10(d.compress(w)), deg=6)
        d = n.where(d == 0, 0, d - 10**n.polyval(poly,cfreqs))
    window = a.dsp.gen_window(d.size, window='kaiser3')
    #window = a.dsp.gen_window(d.size, window='hanning')
    #window = 1
    _d = n.fft.ifft(d * window)
    _w = n.fft.ifft(w * window)
    tau = n.fft.fftfreq(d.size, sdf)
    g = n.sqrt(n.average((w*window)**2))
    _d_c,info = a.deconv.clean(_d, _w, tol=1e-4)
    _d_c += info['res'] / g
    d_c = _d_c.copy()
    #d_c[5:-4] = 0
    wgt = n.exp(-n.abs(n.abs(tau)-1200)**2 / 200**2)
    d_c *= wgt
    #P.subplot(122); P.semilogy(tau, n.abs(d_c), '.')
    #d_c = n.where(n.abs(n.abs(tau)-1200) > 200, 0, d_c)
    d_c = n.fft.fft(d_c)
    j = i+1
    P.figure(1)
    P.subplot(nrows,ncols,j)
    if False:
        P.semilogy(cfreqs, n.abs(d))
        P.semilogy(cfreqs, n.abs(d_c))
        #P.semilogy(cfreqs, n.abs(d_df))
    else:
        P.plot(cfreqs, d, ',')
        P.plot(cfreqs, d_c, ',')
        #P.plot(cfreqs, d_df)
    P.xlabel('GHz')
    P.figure(2)
    P.subplot(nrows,ncols,j)
    P.title(a.miriad.bl2ij(bl))
    #P.bar(tau, n.log10(n.abs(_d_c))+4, width=tau[1]-tau[0])
    P.semilogy(tau, n.abs(_d_c))
    P.xlabel('ns')
    #P.loglog(n.abs(tau), n.abs(_d_c), '.')
P.show()
