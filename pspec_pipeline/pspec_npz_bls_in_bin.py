#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys, optparse

WINDOW = 'blackman-harris'

o = optparse.OptionParser()
o.add_option('-b','--bin', dest='bin', type='int',
    help='The bin to query.')
opts,args = o.parse_args(sys.argv[1:])

dat,wgt = 0, 0
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    bins = f['bins']
    freqs = f['freqs']
    try: cnt = n.argwhere(bins == opts.bin)[0].squeeze()
    except(IndexError):
        print '    bin not found.  Skipping...'
        continue
    bls = f['bls'][cnt]
    u,v,lst = C.pspec.bin2uv(opts.bin)
    print '   ', opts.bin
    print '   ', u, v, lst
    print '    UMAG:', n.sqrt(u**2 + v**2)
    print '   ', [a.miriad.bl2ij(bl) for bl in bls if bl != 0]
    d,w = f['dat'][cnt], f['wgt'][cnt]
    dat += d; wgt += w

for d,w in zip(dat,wgt):
    d = d / n.where(w > 0, w, 1)
    p.subplot(311)
    p.plot(freqs, w)

d, w = dat.sum(axis=0), wgt.sum(axis=0)
window = a.dsp.gen_window(freqs.size, window=WINDOW)
_d = n.fft.ifft(d*window)
_w = n.fft.ifft(w*window)
_d_cl, info = a.deconv.clean(_d, _w, tol=1e-9, stop_if_div=False)
d_mdl = n.fft.fft(_d_cl)
d_res = (d - d_mdl * w) / n.where(w > 0, w, 1)
d /= n.where(w > 0, w, 1)

p.subplot(311)
p.plot(freqs, w, 'k-')

p.subplot(312)
p.plot(freqs, d_mdl.real, 'b-')
p.plot(freqs, d_mdl.imag, 'r-')
p.plot(freqs, d_res.real, 'b-')
p.plot(freqs, d_res.imag, 'r-')

p.subplot(313)
B = .008
kwargs = {'cen_fqs':[.160],'B':B, 'ntaps':3, 'window':WINDOW, 'bm_fqs':freqs.clip(.120,.190)}
Tmdl,ks = C.pspec.Trms_vs_fq(freqs, d_mdl, **kwargs)
for d,w in zip(dat, wgt):
    try: wgt = w.max()
    except: continue
    d = d / wgt; w = w / wgt
    d_res = d - d_mdl * w
    p.subplot(312)
    d_res_wgt = d_res / n.where(w > 0, w, 1)
    p.plot(freqs, d_res_wgt.real, 'b,', alpha=.25)
    p.plot(freqs, d_res_wgt.imag, 'r,', alpha=.25)
    p.subplot(313)
    Tres,ks = C.pspec.Trms_vs_fq(freqs, d_res, **kwargs)
    W,ks = C.pspec.Trms_vs_fq(freqs, w, **kwargs)

    for fq in Tres:
        T_cl, info = a.deconv.clean(Tres[fq], W[fq], tol=1e-9, maxiter=100,
            stop_if_div=False, verbose=False)
        T = Tmdl[fq] + T_cl + info['res']
        p.plot(T, 'k.')
        p.plot(Tmdl[fq], 'g.')
        p.plot(T_cl, 'b.')
        p.plot(info['res'], 'r.')
    print fq

p.show()

