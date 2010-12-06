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
DO_K3PK = True

DECONV1WIDTH = 10
NSIG = 1.5
PLOT2 = False
NDAY = 17

# Read data from npz files (which are already averaged over an hour)
dat,wgt = {},{}
for filename in args:
    print 'Reading', filename
    # Fudge factor to avoid recalibrating
    factor = 1.
    if float('.'.join(filename.split('.')[1:3])) > 2455022.6: factor = 1.5**2
    buf = n.load(filename)
    freqs = buf['freqs']
    for dbl in buf.files:
        if not dbl.startswith('d_'): continue
        #if not a.miriad.bl2ij(int(dbl[2:]))[0] in [0,4,10,12]: continue
        wbl = 'w'+dbl[1:]
        bl = int(dbl[2:])
        if not dat.has_key(bl): dat[bl],wgt[bl] = [],[]
        dat[bl].append(buf[dbl]/factor)
        wgt[bl].append(buf[wbl])

# Put data in temp units and remove smooth models
for bl in dat:
    dat[bl],wgt[bl] = SCALE*n.array(dat[bl]), n.array(wgt[bl])
    buf = []
    for d,w in zip(dat[bl],wgt[bl]):
        # Fudge calibration by using an approximate scaling
        val = n.where(d == 0, 0, 1)
        gain = n.sqrt(n.average(w**2))
        if gain != 0:
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
        buf.append(d*val)
    dat[bl] = n.array(buf)[:,MIN_CH:MAX_CH] #* 2*(n.random.randint(2)-.5)
    wgt[bl] = wgt[bl][:,MIN_CH:MAX_CH]
    # Remove smooth-in-time components
    #avg = n.sum(dat[bl], axis=0) / n.sum(wgt[bl], axis=0).clip(1,n.Inf)
    #dat[bl] -= avg * wgt[bl]
    for ch in range(dat[bl].shape[1]):
        d,w = dat[bl][:,ch], wgt[bl][:,ch]
        gain = n.sqrt(n.average(w**2))
        if gain == 0: continue
        _d = n.fft.ifft(d)
        ker = n.fft.ifft(w)
        __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
        __d[1+NDAY:-NDAY] = 0 # This must match # of days of files provided
        dat[bl][:,ch] -= n.fft.fft(__d) * n.fft.fft(ker)
    if False: # Decorrelate the data for noise estimates
        dat[bl] *= 2*(n.random.randint(2,size=dat[bl].shape)-.5)

# Create data cubes (bl,time,fq) for flagging data statistically along any axis
freqs = freqs[MIN_CH:MAX_CH]
dcube, wcube = [], []
blorder = []
for bl in dat:
    blorder.append(bl)
    dcube.append(dat[bl])
    wcube.append(wgt[bl])
dcube, wcube = n.array(dcube), n.array(wcube)
del(dat); del(wgt)
print dcube.shape
valid = n.ones_like(dcube)
for ax in [0,1,2]:
    _avg = dcube.sum(axis=ax) / wcube.sum(axis=ax).clip(1,n.Inf)
    sh = n.array(dcube.shape); sh[ax] = 1
    dif = dcube - _avg.reshape(sh) * wcube
    sig = n.sqrt((dif**2).sum(axis=ax) / (wcube**2).sum(axis=ax).clip(1,n.Inf))
    sig = sig.reshape(sh)
    valid *= n.where(n.abs(dif) < NSIG*sig*wcube.clip(1,n.Inf), 1, 0)
dcube *= valid
wcube *= valid

# Sum times into different bins
bins = [int(f.split('.')[2][0]) for f in args]
dat,wgt = {},{}
for cnt1,bl in enumerate(blorder):
    dbuf,wbuf = {}, {}
    for cnt2,bin in enumerate(bins):
        if not dat.has_key(bin): dat[bin],wgt[bin] = {},{}
        dat[bin][bl] = dat[bin].get(bl,0) + dcube[cnt1][cnt2]
        wgt[bin][bl] = wgt[bin].get(bl,0) + wcube[cnt1][cnt2]
        #dbuf[bin] = dbuf.get(bin,0) + dcube[cnt1][cnt2]
        #wbuf[bin] = wbuf.get(bin,0) + wcube[cnt1][cnt2]
    # Just examine a single bin for now
    #dat[bl],wgt[bl] = dbuf[4],wbuf[4]

## Plot data
#nplts = 3
#m1 = n.ceil(n.sqrt(nplts))
#m2 = n.ceil(nplts / m1)
#if True:
#    nplts += len(blorder)
#    m1 = n.ceil(n.sqrt(nplts))
#    m2 = n.ceil(nplts / m1)
#    for i, bl in enumerate(blorder):
#        d,w = dcube[i], wcube[i]
#        if True:
#            p.subplot(m2, m1, i+1)
#            i,j = a.miriad.bl2ij(bl)
#            #_d,_w = n.fft.rfft(d,axis=0), n.fft.rfft(w,axis=0)
#            #p.imshow(n.log10(n.abs(_d) / n.abs(_w[0:1,:]).clip(1,n.Inf)), vmax=3, vmin=1, aspect='auto')
#            p.imshow(n.log10(n.abs(d / w.clip(1,n.Inf))), vmax=4, vmin=2, aspect='auto')
#            p.title(str(i))
#    p.show()

# Compute the median per-antenna score and flag off antennas above it
if False:
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
__dat,__wgt = {},{}
for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
    ch0 = b*SUBBAND/2
    ch1 = ch0 + NTAPS*SUBBAND
    if ch1 > MAX_CH-MIN_CH: continue
    for lstbin in dat:
        if not __dat.has_key(lstbin): __dat[lstbin],__wgt[lstbin] = {},{}
        if not data.has_key(lstbin): data[lstbin],gain[lstbin] = {},{}
        for i in dat[lstbin]:
            _d = C.pfb.pfb(dat[lstbin][i][ch0:ch1], taps=NTAPS, window='kaiser3', fft=n.fft.ifft)
            _w = C.pfb.pfb(wgt[lstbin][i][ch0:ch1], taps=NTAPS, window='kaiser3', fft=n.fft.ifft)
            g1 = n.abs(_w[0])
            __d, info = a.deconv.clean(_d, _w, tol=opts.clean)
            __d[1:] = 0
            _d -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(_w))
            _d1 = _d[:_d.size/2]
            # Cross correlate with other PS measurements
            for _d2,g2 in zip(data[lstbin].get(b,[]), gain[lstbin].get(b,[])):
                __dat[lstbin][b] = __dat[lstbin].get(b,0) + _d1 * n.conj(_d2)
                # Note that as of here wgt is just a number, not a vector
                __wgt[lstbin][b] = __wgt[lstbin].get(b,0) + g1 * g2
            data[lstbin][b] = data[lstbin].get(b,[]) + [_d1]
            gain[lstbin][b] = gain[lstbin].get(b,[]) + [g1]

# Convert into cosmological units (Mpc, mK^2, etc)
dat = {}
dtot,wtot = 0,0
for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
    ch0 = b*SUBBAND/2
    ch1 = ch0 + NTAPS*SUBBAND
    if ch1 > MAX_CH-MIN_CH: continue
    _freqs = freqs[ch0+int((NTAPS-1)/2.):ch0+int((NTAPS-1)/2.)+SUBBAND]
    B = (_freqs[-1] - _freqs[0])
    z = C.pspec.f2z(_freqs)
    eta = C.pspec.f2eta(_freqs)
    k_pl = (eta * C.pspec.dk_deta(z))[:int(n.ceil(float(eta.size)/2))]
    k_pl = k_pl.clip(k_pl[1]/2,n.Inf) # Having k=0 messes up log binning
    if DO_K3PK: scalar = C.pspec.k3pk_from_Trms(1, k_pl, n.average(_freqs), B=B)
    else: scalar = 1
    for lstbin in __dat:
        if not dat.has_key(lstbin): dat[lstbin] = []
        ks, d = C.pspec.rebin_log(k_pl, __dat[lstbin][b] * scalar, nbins=8)
        dtot += d
        wtot += __wgt[lstbin][b]
        dat[lstbin].append(d / __wgt[lstbin][b])

# Plot the data
if DO_K3PK: mx,drng = 11,6
else: mx,drng = 4,4
#mx,drng = dat1.max(), dat1.max()-max(dat1.min(),2)
nplts = len(dat) + 1
m1 = int(n.ceil(n.sqrt(nplts)))
m2 = int(n.ceil(float(nplts)/m1))
for cnt, lstbin in enumerate(dat):
    dat[lstbin] = n.array(dat[lstbin]).real
    p.subplot(m2, m1, cnt+1)
    d = n.log10(n.abs(dat[lstbin]))
    p.imshow(d, vmax=mx, vmin=mx-drng, aspect='auto', interpolation='nearest')
p.subplot(m2,m1, cnt+2)
dtot = n.abs(dtot.real)
print ks
print dtot/wtot
p.loglog(ks, n.abs(dtot/wtot).clip(1,n.Inf))
p.grid()
p.show()


