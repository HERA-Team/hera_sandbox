#!/usr/bin/python
import aipy as a, numpy as n, pylab as p, capo as C
import sys, optparse, os

o = optparse.OptionParser()
o.set_usage(os.path.basename(sys.argv[0]) + ' [options] *autoex.npz')
o.set_description(__doc__)
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

MIN_CH,MAX_CH,SUBBAND = 250,750,50 # For 1024ch GB data

#DECONV1WIDTH = 5
DECONV1WIDTH = 10
#NSIG = 2
NSIG = 1.5
#NSIG = 1
PLOT1 = False
PLOT2 = False

# Read in data
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
        if not a.miriad.bl2ij(int(dbl[2:]))[0] in [0,4,10,12]: continue
        wbl = 'w'+dbl[1:]
        bl = int(dbl[2:])
        if not dat.has_key(bl): dat[bl],wgt[bl] = [],[]
        dat[bl].append(buf[dbl]/factor)
        wgt[bl].append(buf[wbl])

# Put data in temp units and remove smooth models
for bl in dat:
    print a.miriad.bl2ij(bl)[0]
    dat[bl],wgt[bl] = n.array(dat[bl]), n.array(wgt[bl])
    # Fudge calibration by assuming 500K median Tsys across band
    avg = dat[bl] / wgt[bl].clip(1,n.Inf)
    dat[bl] *= (500 / n.median(avg))
    wgt[bl] *= n.where(dat[bl] == 0, 0, 1)
    buf = []
    for d,w in zip(dat[bl],wgt[bl]):
        val = n.where(d == 0, 0, 1)
        gain = n.sqrt(n.average(w**2))
        if gain != 0:
            # 1st pass: deconv by power law model across entire band
            poly = n.polyfit(n.log10(freqs.compress(val)), n.log10(n.abs(d.compress(val))), deg=1)
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
    dat[bl] = n.array(buf)[:,MIN_CH:MAX_CH]
    wgt[bl] = wgt[bl][:,MIN_CH:MAX_CH]

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

# Plot data
nplts = 3
m1 = n.ceil(n.sqrt(nplts))
m2 = n.ceil(nplts / m1)
if False:
    nplts += len(blorder)
    m1 = n.ceil(n.sqrt(nplts))
    m2 = n.ceil(nplts / m1)
    for i, bl in enumerate(blorder):
        d,w = dcube[i], wcube[i]
        if True:
            p.subplot(m2, m1, i+1)
            i,j = a.miriad.bl2ij(bl)
            p.imshow(n.log10(n.abs(d / w.clip(1,n.Inf))), vmax=0.5, vmin=-1.5, aspect='auto')
            p.title(str(i))
if True:
    p.subplot(m2,m1,nplts-2)
    ds,ws = dcube.sum(axis=0), wcube.sum(axis=0)
    avg = ds / ws.clip(1,n.Inf)
    if False:
        _avg = n.fft.ifft2(avg)
        _avg = n.concatenate([_avg[_avg.shape[0]/2:], _avg[:_avg.shape[0]/2]], axis=0)
        _avg = n.concatenate([_avg[:,_avg.shape[1]/2:], _avg[:,:_avg.shape[1]/2]], axis=1)
        p.imshow(n.log10(n.abs(_avg)), vmax=0, vmin=-3, aspect='auto')
    else:
        p.imshow(n.log10(n.abs(ds) / ws.clip(1,n.Inf)), vmax=0.5, vmin=-1.5, aspect='auto')
    p.colorbar()

    p.subplot(m2,m1,nplts-1)
    dss,wss = ds.sum(axis=0), ws.sum(axis=0)
    spec = dss / wss.clip(1,n.Inf)
    #p.plot(spec)

    p.subplot(m2,m1,nplts)
    _spec = n.fft.ifft(spec)
    p.semilogy(n.abs(n.concatenate([_spec[_spec.size/2:], _spec[:_spec.size/2]])))

    p.subplot(m2,m1,nplts-1)
    _d = n.fft.ifft(dss)
    ker = n.fft.ifft(wss)
    __d, info = a.deconv.clean(_d, ker, tol=opts.clean)
    a__d = n.abs(__d)
    __d = n.where(a__d > 0.25 * a__d.max(), __d, 0)
    mdl = n.fft.fft(__d) * wss
    dss -= mdl.real
    spec = dss / wss.clip(1,n.Inf)
    p.plot(spec)
    p.ylim(-0.5,0.5)

    p.subplot(m2,m1,nplts)
    p.semilogy(n.abs(n.concatenate([__d[__d.size/2:], __d[:__d.size/2]])), '^')
    _spec = n.fft.ifft(spec)
    p.semilogy(n.abs(n.concatenate([_spec[_spec.size/2:], _spec[:_spec.size/2]])))

p.show()

