#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys, optparse, os

o = optparse.OptionParser()
o.set_usage(os.path.basename(sys.argv[0]) + ' [options] *autoex.npz')
o.set_description(__doc__)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

MIN_CH,MAX_CH,SUBBAND = 250,750,50  # for 1024ch GB data
NTAPS = 3

#SCALE = 30e3
SCALE = 6*30e3
DO_K3PK = False

#DECONV1WIDTH = 5
DECONV1WIDTH = 10
#NSIG = 2
NSIG = 1.5
#NSIG = 1
PLOT2 = False

# Read data from npz files (which are already averaged over an hour)
dat,wgt = {},{}
for filename in args:
    print 'Reading', filename
    buf = n.load(filename)
    freqs = buf['freqs']
    for dbl in buf.files:
        if not dbl.startswith('d_'): continue
        wbl = 'w'+dbl[1:]
        # Fudge calibration by using an approximate scaling
        d = buf[dbl] * SCALE
        w = buf[wbl]
        val = n.where(d == 0, 0, 1)
        gain = n.sqrt(n.average(w**2))
        if gain == 0: continue
        # 1st pass: deconv by power law model across entire band
        poly = n.polyfit(n.log10(freqs.compress(val)), n.log10(n.abs(d.compress(val))), deg=3)
        mdl = 10**(n.polyval(poly, n.log10(freqs)))
        ker = n.fft.ifft(w*mdl)
        _d = n.fft.ifft(d)
        __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
        if DECONV1WIDTH == 0: __d[1:] = 0
        else: __d[DECONV1WIDTH+1:-DECONV1WIDTH] = 0
        d -= n.fft.fft(__d) * n.fft.fft(ker)
        # 2nd pass: deconv across entire band (no power law)
        _d = n.fft.ifft(d)
        ker = n.fft.ifft(w)
        __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
        if DECONV1WIDTH == 0: __d[1:] = 0
        else: __d[DECONV1WIDTH+1:-DECONV1WIDTH] = 0
        d -= n.fft.fft(__d) * n.fft.fft(ker)
        # 3rd pass: deconv only in band of interest
        _d = n.fft.ifft(d[MIN_CH:MAX_CH])
        ker = n.fft.ifft(w[MIN_CH:MAX_CH])
        __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
        if DECONV1WIDTH == 0: __d[1:] = 0
        else: __d[DECONV1WIDTH+1:-DECONV1WIDTH] = 0
        d[MIN_CH:MAX_CH] -= n.fft.fft(__d) * n.fft.fft(ker)
        i = a.miriad.bl2ij(int(dbl[2:]))[0]
        #dat[i] = dat.get(i,0) + d[MIN_CH:MAX_CH]
        dat[i] = dat.get(i,0) + d[MIN_CH:MAX_CH] * 2*(n.random.randint(2)-.5)
        wgt[i] = wgt.get(i,0) + w[MIN_CH:MAX_CH]
freqs = freqs[MIN_CH:MAX_CH]

# Compute the median per-antenna score and flag off antennas above it
scores = {}
for i in dat:
    avg = dat[i] / wgt[i].clip(1,n.Inf)
    scores[i] = n.std(avg[100:-100])
mscore = n.median(scores.values())
for i in scores:
    print i, scores[i],
    if scores[i] > mscore:
        print '*'
        del(dat[i])
        del(wgt[i])
        continue
    print ''

if PLOT2:
    totavg = sum([dat[i] for i in dat]) / sum([wgt[i] for i in wgt]).clip(1,n.Inf)
    p.subplot(121)
    p.plot(freqs, totavg, label='avg')
    p.ylim(-10e3,10e3)
    p.subplot(122)
    pspec = n.fft.ifft(totavg)
    pspec = n.concatenate([pspec[pspec.size/2:], pspec[:pspec.size/2]])
    p.semilogy(n.abs(pspec), label='avg')
    p.legend()
    p.show()

# Perform line-of-sight transform on each subband and cross-multiply between antennas
data,gain = {},{}
__dat,__wgt,__dec,__cln = {},{},{},{}
for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
    ch0 = b*SUBBAND/2
    ch1 = ch0 + NTAPS*SUBBAND
    if ch1 > MAX_CH-MIN_CH: continue
    for i in dat:
        _d = C.pfb.pfb(dat[i][ch0:ch1], taps=NTAPS, window='kaiser3', fft=n.fft.ifft)
        _w = C.pfb.pfb(wgt[i][ch0:ch1], taps=NTAPS, window='kaiser3', fft=n.fft.ifft)
        g1 = n.abs(_w[0])
        __d, info = a.deconv.clean(_d, _w, tol=opts.clean)
        __d[1:] = 0
        _d -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(_w))
        _d1 = _d[:_d.size/2]
        # Cross correlate with other PS measurements
        for _d2,g2 in zip(data.get(b,[]), gain.get(b,[])):
            __dat[b] = __dat.get(b,0) + _d1 * n.conj(_d2)
            __wgt[b] = __wgt.get(b,0) + g1 * g2
        data[b] = data.get(b,[]) + [_d1]
        gain[b] = gain.get(b,[]) + [g1]

# Convert into cosmological units (Mpc, mK^2, etc)
dat = []
for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
    ch0 = b*SUBBAND/2
    ch1 = ch0 + NTAPS*SUBBAND
    if ch1 > MAX_CH-MIN_CH: continue
    _freqs = freqs[ch0+int((NTAPS-1)/2.):ch0+int((NTAPS-1)/2.)+SUBBAND]
    B = (_freqs[-1] - _freqs[0])
    z = C.pspec.f2z(_freqs)
    eta = C.pspec.f2eta(_freqs)
    k_pl = (eta * C.pspec.dk_deta(z))[:__dat[b].size]
    if DO_K3PK: scalar = C.pspec.k3pk_from_Trms(1, k_pl, n.average(_freqs), B=B)
    else: scalar = 1
    dat.append((__dat[b] / __wgt[b]) * scalar)
dat = n.array(dat).real
#avg = n.average(dat,axis=0); avg.shape = (1,avg.size); dat -= avg

#p.subplot(221)
dat1 = n.log10(n.abs(dat))
#mx,drng = 10,7
#mx,drng = 7,5
#mx,drng = 3,4
if DO_K3PK: mx,drng = 10,4
else: mx,drng = 6,6
#mx,drng = dat1.max(), dat1.max()-max(dat1.min(),2)
p.imshow(dat1, vmax=mx, vmin=mx-drng, aspect='auto', interpolation='nearest')
p.colorbar(shrink=.5)
p.show()


