#! /usr/bin/env python
import aipy as a, numpy as n, capo as C
import optparse, sys
import pylab as p

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True)
o.add_option('--sep', dest='sep',
    help='Separation to use.')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
opts,args = o.parse_args(sys.argv[1:])

freqs = n.linspace(.1, .2, 203)
chans = a.scripting.parse_chans(opts.chan, freqs.size)
WINDOW = opts.window
sep = opts.sep

afreqs = freqs.take(chans)
fq = n.average(afreqs)
z = C.pspec.f2z(fq)
sdf = freqs[1]-freqs[0]

t = n.arange(-200,200,1) * 42.8
w = a.dsp.gen_window(t.size, 'none')
sig = .000288*(fq/.1788)
cen = .001059 * (fq/.1788)
fir = n.exp(-(t**2)/(2*(1./sig)**2)).astype(n.complex) * w
fir /= fir.sum()
fir *= n.exp(2j*n.pi*cen*t) # need to flip the sign for backward baselines
fir = n.array([1.])
p.plot(fir.real)
p.plot(fir.imag)
p.show()

B = sdf * afreqs.size * C.pfb.NOISE_EQUIV_BW[WINDOW]
etas = n.fft.fftshift(C.pspec.f2eta(afreqs))
kpl = etas * C.pspec.dk_deta(z)
bm = n.polyval(C.pspec.DEFAULT_BEAM_POLY, fq)
scalar = C.pspec.X2Y(z) * bm * B

print 'Freq:',fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar
print 'sep:', sep

_T, _W = {}, {}
for filename in args:
    print 'Reading', filename
    npz = n.load(filename)
    d = npz[sep].take(chans, axis=1)
    try: w = npz['wgt_'+sep].take(chans, axis=1) / 1e-6
    except(KeyError): w = n.ones_like(d)
    d = n.concatenate([d[-200:],d[:1460]], axis=0)
    w = n.concatenate([w[-200:],w[:1460]], axis=0)
    #w = n.where(w > .1, 1., 0)
    for ch in xrange(d.shape[1]):
        d[:,ch] = n.convolve(w[:,ch]*d[:,ch], fir, mode='same')
        w[:,ch] = n.convolve(w[:,ch], n.abs(fir), mode='same')
        d[:,ch] /= n.where(w[:,ch] > 0, w[:,ch], 1)
    #w = n.where(w > .1, 1, 0)
    #w = n.average(w, axis=1).reshape((w.shape[0],1)) * n.ones((1,w.shape[1]))
    p.subplot(131); C.arp.waterfall(d, mode='phs')
    p.subplot(132); C.arp.waterfall(w)
    d *= w
    #C.arp.waterfall(n.fft.fft(d[300:1000], axis=0), mode='lin')
    #p.colorbar()
    #p.show()
    T = d * C.pspec.jy2T(afreqs)
    window = a.dsp.gen_window(T.shape[1], WINDOW); window.shape = (1,window.size)
    _T[filename] = n.fft.fftshift(n.fft.ifft(window * T), axes=1)
    #_W[filename] = n.fft.fftshift(n.fft.ifft(window * w), axes=1)
    _W[filename] = n.fft.ifft(window * w)

    p.subplot(133); C.arp.waterfall(_T[filename], drng=3)
    #p.subplot(122); C.arp.waterfall(_W[filename], drng=3)
    p.show()

import pylab as p

for i,filename in enumerate(_T.keys()):
    print i, filename
    p.subplot(1, len(_T), i+1)
    C.arp.waterfall(_T[filename], mode='log')
    p.colorbar(shrink=.5)
p.show()

pk_sum, pk_wgt = 0., 0.
for i,f1 in enumerate(_T.keys()):
    _T1 = _T[f1]
    _W1 = _W[f1]
    for f2 in _T.keys()[i+1:]:
        _T2 = _T[f2]
        _W2 = _W[f2]
        pk12 = scalar * _T1 * n.conj(_T2)
        wt12 = _W1 * _W2
        wt12 /= wt12.max()
        pk_sum += pk12
        p.subplot(122)
        p.plot(kpl, n.average(pk12[200:1200], axis=0).real)
        #pk_wgt += 1
        #p.clf(); p.plot(wt12[:,0]); p.show()
        pk_wgt += wt12[:,:1]

p.subplot(121)
pk = pk_sum / pk_wgt
#C.arp.waterfall(pk, mode='real', mx=1e8, drng=2e8)#, drng=3)
C.arp.waterfall(pk, mode='real', mx=1e11, drng=2e11)#, drng=3)
p.colorbar()

pk_avg = pk_sum[200:1200].sum(axis=0) / pk_wgt[200:1200].sum(axis=0)
print pk_avg.shape

p.subplot(122)
p.plot(kpl, pk_avg.real, 'k^-')
p.show()

p.plot(n.abs(kpl), 2* n.abs(kpl)**3/(2*n.pi**2) * pk_avg.real, 'k.-')
p.show()

_T = n.concatenate([_T[f] for f in _T.keys()], axis=1)

_T = _T[200:1200].T
C.arp.waterfall(cov(_T), mode='log')
p.show()


