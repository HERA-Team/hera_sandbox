#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys

def apply_lpf(d, w, width, iter, lpf=None):
    N = d.size
    if lpf is None:
        _lpf = n.ones(N, dtype=n.float)
        _lpf[width+1:-width] = 0
        lpf = n.fft.fftshift(n.fft.fft(_lpf))
    if iter == 0:
        d_lpf = n.convolve(d, lpf, 'full')[N/2:3*N/2] / N
    else:
        if True:
            d_lpf = 0
            for i in range(iter):
                d_lpf += apply_lpf((d-d_lpf)*w, w, width, 0, lpf=lpf)
        else:
            d_lpf = apply_lpf(d, w, width, iter-1, lpf=lpf)
            w_lpf = apply_lpf(d_lpf*w, w, width, iter-1, lpf=lpf)
            correction = d_lpf / w_lpf
            d_lpf *= correction
    return d_lpf


for f in sys.argv[1:]:
    uv = a.miriad.UV(f)
    cnt = 0
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        cnt += 1
        if cnt < 41: continue
        window = a.dsp.gen_window(d.size, 'blackman-harris')
        w = n.logical_not(f).astype(n.float)
        d *= w
        d_lpf = apply_lpf(d,w,18,10)
        d_rs = (d - d_lpf) * w
        p.subplot(131); p.plot(d_lpf*w)
        p.subplot(132); p.semilogy(n.fft.fftshift(n.abs(n.fft.ifft(d_rs))))
        p.subplot(133); p.semilogy(n.fft.fftshift(n.abs(n.fft.ifft(d_rs*window))))

        for d_ in [d,d_rs]:
            _d,_w = n.fft.ifft(d_*window), n.fft.ifft(w*window)
            _d,info = a.deconv.clean(_d, _w, tol=1e-5)
            d_cl = n.fft.fft(_d)
            d_rs = (d_ - d_cl) * w
            p.subplot(131); p.plot(d_cl*w)
            p.subplot(132); p.semilogy(n.fft.fftshift(n.abs(n.fft.ifft(d_rs))))
            p.subplot(133); p.semilogy(n.fft.fftshift(n.abs(_d+info['res'])))

        #p.subplot(131); p.plot(d)
        #p.subplot(132); p.semilogy(n.fft.fftshift(n.abs(n.fft.ifft(d))))
        #p.subplot(133); p.semilogy(n.fft.fftshift(n.abs(n.fft.ifft(d*window))))
        #p.subplot(224)
        #p.semilogy(n.fft.fftshift(n.abs(n.fft.ifft(d_rs))))
            
        p.show()
