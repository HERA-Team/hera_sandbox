#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys, optparse, os

o = optparse.OptionParser()
o.set_usage(os.path.basename(sys.argv[0]) + ' [options] *cross.npz')
o.set_description(__doc__)
o.add_option('--pk', dest='pk', action='store_true',
    help='Plot P(k) instead of k^3 P(k)')
o.add_option('--cutoff', dest='cutoff', type='float', default=1,
    help='Number of integerations required to be in a bin that is processed.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

MIN_CH,MAX_CH,SUBBAND = 250,700,50  # for 1024ch GB data
NTAPS = 3

DECONV1WIDTH = 10
NSIG = 1.5
PLOT2 = False

# Read data from npz files (which are already averaged over an hour)
dat,wgt = {},{}
for filename in args:
    print 'Reading', filename
    buf = n.load(filename)
    freqs = buf['freqs']
    times = buf['times']
    for dbl in buf.files:
        if not dbl.startswith('d_'): continue
        wbl = 'w'+dbl[1:]
        bin = int(dbl[2:])
        if not dat.has_key(bin): dat[bin],wgt[bin] = [],[]
        dat[bin].append(buf[dbl])
        wgt[bin].append(buf[wbl])

SCALE = C.pspec.jy2T(freqs)
# Put data in temp units and remove smooth models
datbins = dat.keys()
for bin in datbins:
    dat[bin],wgt[bin] = SCALE*n.array(dat[bin]), n.array(wgt[bin])
    m = n.max(wgt[bin].sum(axis=0))
    if m < opts.cutoff:
        del(dat[bin]); del(wgt[bin])
        continue
    else:
        print m,
        print bin, C.pspec.bin2uv(bin)
    print n.abs(dat[bin][:,MIN_CH:MAX_CH]).sum()/wgt[bin][:,MIN_CH:MAX_CH].sum()
    buf = []
    for d,w in zip(dat[bin],wgt[bin]):
        # Fudge calibration by using an approximate scaling
        val = n.where(d == 0, 0, 1)
        gain = n.sqrt(n.average(w**2))
        if gain != 0:
        #if False:
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
    dat[bin] = n.array(buf)[:,MIN_CH:MAX_CH] #* 2*(n.random.randint(2)-.5)
    wgt[bin] = wgt[bin][:,MIN_CH:MAX_CH]
    print n.abs(dat[bin]).sum()/wgt[bin].sum()
    if False:
        d = dat[bin].sum(axis=0)
        w = wgt[bin].sum(axis=0).clip(1,n.Inf)
        pltdat = n.abs(d/w).squeeze()
        p.semilogy(pltdat, '.')
        #p.loglog(w, n.abs(d/w).clip(.001,n.Inf), '.')
    ## Remove smooth-in-time components
    ##avg = n.sum(dat[bin], axis=0) / n.sum(wgt[bin], axis=0).clip(1,n.Inf)
    ##dat[bin] -= avg * wgt[bin]
    #for ch in range(dat[bin].shape[1]):
    #    d,w = dat[bin][:,ch], wgt[bin][:,ch]
    #    gain = n.sqrt(n.average(w**2))
    #    if gain == 0: continue
    #    _d = n.fft.ifft(d)
    #    ker = n.fft.ifft(w)
    #    __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
    #    __d[1+NDAY:-NDAY] = 0 # This must match # of days of files provided
    #    dat[bin][:,ch] -= n.fft.fft(__d) * n.fft.fft(ker)
    if False: # Decorrelate the data for noise estimates
        dat[bin] *= 2*(n.random.randint(2,size=dat[bin].shape)-.5)

#p.show()

## Create data cubes (bin,time,fq) for flagging data statistically along any axis
#freqs = freqs[MIN_CH:MAX_CH]
#dcube, wcube = [], []
#binorder = []
#for bin in dat:
#    binorder.append(bin)
#    dcube.append(dat[bin])
#    wcube.append(wgt[bin])
#dcube, wcube = n.array(dcube), n.array(wcube)
#del(dat); del(wgt)
#print dcube.shape
#valid = n.ones_like(dcube)
#for ax in [0,1,2]:
#    _avg = dcube.sum(axis=ax) / wcube.sum(axis=ax).clip(1,n.Inf)
#    sh = n.array(dcube.shape); sh[ax] = 1
#    dif = dcube - _avg.reshape(sh) * wcube
#    sig = n.sqrt((dif**2).sum(axis=ax) / (wcube**2).sum(axis=ax).clip(1,n.Inf))
#    sig = sig.reshape(sh)
#    valid *= n.where(n.abs(dif) < NSIG*sig*wcube.clip(1,n.Inf), 1, 0)
#dcube *= valid
#wcube *= valid

## Plot data
#nplts = 3
#m1 = n.ceil(n.sqrt(nplts))
#m2 = n.ceil(nplts / m1)
#if True:
#    nplts += len(binorder)
#    m1 = n.ceil(n.sqrt(nplts))
#    m2 = n.ceil(nplts / m1)
#    for i, bin in enumerate(binorder):
#        d,w = dcube[i], wcube[i]
#        if True:
#            p.subplot(m2, m1, i+1)
#            i,j = a.miriad.bin2ij(bin)
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
    for bin in dat:
        if not __dat.has_key(bin): __dat[bin],__wgt[bin] = {},{}
        if not data.has_key(bin): data[bin],gain[bin] = {},{}
        for i in range(dat[bin].shape[0]):
            _d = C.pfb.pfb(dat[bin][i,ch0:ch1], taps=NTAPS, window='kaiser3', fft=n.fft.ifft)
            _w = C.pfb.pfb(wgt[bin][i,ch0:ch1], taps=NTAPS, window='kaiser3', fft=n.fft.ifft)
            g1 = n.abs(_w[0])
            if g1 != 0:
                __d, info = a.deconv.clean(_d, _w, tol=opts.clean)
                __d[1:] = 0
                _d -= n.fft.ifft(n.fft.fft(__d) * n.fft.fft(_w))
                _d1 = _d[:_d.size/2]
            else: _d1 = n.zeros(_d.size/2)
            # Cross correlate with other PS measurements
            for _d2,g2 in zip(data[bin].get(b,[]), gain[bin].get(b,[])):
                __dat[bin][b] = __dat[bin].get(b,0) + _d1 * n.conj(_d2)
                # Note that as of here wgt is just a number, not a vector
                __wgt[bin][b] = __wgt[bin].get(b,0) + g1 * g2
            data[bin][b] = data[bin].get(b,[]) + [_d1]
            gain[bin][b] = gain[bin].get(b,[]) + [g1]

# Convert into cosmological units (Mpc, mK^2, etc)
dat = {}
dtot,wtot = [],[]
for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
    ch0 = b*SUBBAND/2
    ch1 = ch0 + NTAPS*SUBBAND
    if ch1 > MAX_CH-MIN_CH: continue
    dtot.append(0); wtot.append(0)
    _freqs = freqs[ch0+int((NTAPS-1)/2.):ch0+int((NTAPS-1)/2.)+SUBBAND]
    B = (_freqs[-1] - _freqs[0])
    z = C.pspec.f2z(_freqs)
    eta = C.pspec.f2eta(_freqs)
    k_pl = (eta * C.pspec.dk_deta(z))[:int(n.ceil(float(eta.size)/2))]
    k_pl = k_pl.clip(k_pl[1]/2,n.Inf) # Having k=0 messes up log binning
    if opts.pk: scalar = 1
    else: scalar = C.pspec.k3pk_from_Trms(1, k_pl, n.average(_freqs), B=B)
    for bin in __dat:
        if not dat.has_key(bin): dat[bin] = []
        ks, d = C.pspec.rebin_log(k_pl, __dat[bin][b] * scalar, nbins=8)
        dtot[-1] += d * __wgt[bin][b]
        wtot[-1] += __wgt[bin][b]**2
        dat[bin].append(d / __wgt[bin][b])
dtot = n.array(dtot)
wtot = n.array(wtot)
wtot.shape = wtot.shape + (1,)

# Plot the data
if opts.pk: mx,drng = 3,4
else: mx,drng = 10,6
#mx,drng = dat1.max(), dat1.max()-max(dat1.min(),2)
nplts = len(dat) + 2
m1 = int(n.ceil(n.sqrt(nplts)))
m2 = int(n.ceil(float(nplts)/m1))
for cnt, bin in enumerate(dat):
    dat[bin] = n.array(dat[bin]).real
    p.subplot(m2, m1, cnt+1)
    d = n.log10(n.abs(dat[bin].real))
    p.imshow(d, vmax=mx, vmin=mx-drng, aspect='auto', interpolation='nearest')
    u,v,lst = C.pspec.bin2uv(bin)
    p.title(str((u,v,'%4.2f'%lst)))

# Plot avg of all bins
p.subplot(m2,m1, cnt+2)
d = n.log10(n.abs(dtot.real/wtot))
p.imshow(d, vmax=mx, vmin=mx-drng, aspect='auto', interpolation='nearest')
p.title('Total')

# Plot avg of all bands
p.subplot(m2,m1, cnt+3)
dtot = n.abs(dtot.real)
for d,w in zip(dtot,wtot):
    p.loglog(ks, n.abs(d/w))
p.loglog(ks, n.abs(dtot.real.sum(axis=0)/wtot.sum()), ':')
p.grid()
p.show()


