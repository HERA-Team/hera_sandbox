#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, math, glob
PLOT = True
if PLOT: import pylab as P; P.ion()

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, chan=True, src=True,
    dec=True, ant=True, pol=True)
o.add_option('-b', '--blsrcs', dest='blsrcs', default='',
    help='Sources (from any tier) that need per-baseline solutions')
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-r', '--drw', dest='drw', type=int, default=5,
    help='The number of delay-rate bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('--maxiter', dest='maxiter', type='int', default=100,
    help='Maximum number of iterations to allow.')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
NCH = chans.size
aa.select_chans(chans)
afreqs = aa.get_afreqs()

srctier = opts.src.split('/')
for i, s in enumerate(srctier):
    srclist, cutoff, catalogs = a.scripting.parse_srcs(s, opts.cat)
    srctier[i] = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
cat = a.fit.SrcCatalog()
for c in srctier: cat.add_srcs(c.values())
blsrcs = opts.blsrcs.split(',')

BURN_IN = 30
CLEAN_GAIN = .2
FRAC_NOISE = 1e-3
NOISE_CLIP = 1e-2
EXP_NOISE = 1.
TRIM = True

# Generate a list of files for padding
def globize(s):
    if s.isdigit(): return '*'
    return s

def gen_filter(shape, dw, drw, ratio=.25):
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
    return filter.clip(0,1)

def coeff_mask(filter, red_dw, red_drw):
    filter_take = n.where(filter)
    red_filter = gen_filter(filter.shape, red_dw, red_drw)
    red_filter_take = n.where(red_filter)
    mask = n.ones_like(filter_take[0])
    for i,(x,y) in enumerate(zip(*filter_take)):
        flag = False
        for _x,_y in zip(*red_filter_take):
            if x == _x and y == _y:
                flag = True
                break
        if not flag: mask[i] = 0
    return mask

def gen_window(shape):
    w = 4 * (.5 - n.abs(.5 - n.arange(shape[0], dtype=n.float)/shape[0])).clip(0,.25)
    w.shape = (w.size, 1)
    return w

try:
    import fftw3
    print 'Using FFTW FFT'
    _fftw3_dat, _fftw3_fwd, _fftw3_rev = None, None, None
    def fft2(d):
        global _fftw3_fwd, _fftw3_dat
        if _fftw3_fwd is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_fwd = fftw3.Plan(_fftw3_dat, None, direction='forward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_fwd()
        return _fftw3_dat
    def ifft2(d):
        global _fftw3_rev, _fftw3_dat
        if _fftw3_rev is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_rev = fftw3.Plan(_fftw3_dat, None, direction='backward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_rev()
        return _fftw3_dat
except(ImportError): 
    print 'Using numpy FFT'
    fft2, ifft2 = n.fft.fft2, n.fft.ifft2

def compute_ant(step):
    if step.has_key('ant'): return
    step['ant'] = {}
    for k in step['srcs']:
        step['ant'][k] = {}
        for i in step['srcs'][k]:
            step['ant'][k][i] = from_coeffs(step['srcs'][k][i])

def compute_residual(msrdat, msrval, simdat, rsvdat, step):
    # Compute residual for all sources
    _resdat = {}
    compute_ant(step)
    for bl in msrdat:
        i,j = a.miriad.bl2ij(bl)
        _resdat[bl] = msrdat[bl].copy()
        for k in simdat:
            if rsvdat.has_key(k): mdl = simdat[k][bl] * rsvdat[k][bl]
            else: mdl = simdat[k][bl]
            try: _resdat[bl] -= step['ant'][k][i] * n.conj(step['ant'][k][j]) * mdl
            except(KeyError): pass
        _resdat[bl] *= msrval[bl] # Mask out invalid data
    step['xtalk'] = compute_xtalk(msrval, _resdat)
    return _resdat

def compute_xtalk(msrval, resdat):
    _xtalk = {}
    for bl in msrval:
        xsum = n.sum(resdat[bl], axis=0)
        xwgt = n.sum(msrval[bl], axis=0)
        _xtalk[bl] = xsum/xwgt.clip(1,n.Inf)
        resdat[bl] -= _xtalk[bl]
        resdat[bl] *= msrval[bl] # Mask out invalid data
    return resdat, _xtalk

file_glob = '.'.join(map(globize, args[0].split('.')))
filelist = glob.glob(file_glob)
filelist.sort()

LINE1 = None
LINE = {}
LINE2 = {}

for arg in args:
    # Gather data
    print 'Processing', arg
    findex = filelist.index(arg)
    files = filelist[findex:findex+1]
    #files = filelist[findex-1:findex+2]
    msrdat, msrval, simdat, rsvdat = {}, {}, {}, {}
    blwgt = {}
    ants = {}
    for k in cat:
        if k in blsrcs: rsvdat[k] = {}
        simdat[k] = {}
        blwgt[k] = {}
    times = []
    # Collect data
    for filename in files:
        print '    Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        uv.select('decimate', opts.decimate, opts.decphs)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            if len(times) == 0 or t != times[-1]:
                times.append(t)
                aa.set_jultime(t)
                cat.compute(aa)
            bl = a.miriad.ij2bl(i,j)
            f = n.logical_not(f.take(chans)).astype(n.int)
            d = d.take(chans) * f / aa.passband(i,j)
            msrdat[bl] = msrdat.get(bl,[]) + [d]
            msrval[bl] = msrval.get(bl,[]) + [f]
            ants[i] = ants[j] = None
            for k in simdat.keys():
                try:
                    simd = n.conj(aa.gen_phs(cat[k],i,j))
                    _f = f
                    u,v,w = aa.gen_uvw(i,j,cat[k]).squeeze()
                    if k in blsrcs:
                        rsvd = aa.resolve_src(u, v, srcshape=cat[k].srcshape)
                        u = v = w = n.ones_like(d)
                except(a.phs.PointingError):
                    _f = n.zeros_like(f)
                    simd = n.zeros_like(d)
                    u = v = w = n.zeros_like(d)
                    if k in blsrcs: rsvd = simd
                simdat[k][bl] = simdat[k].get(bl,[]) + [simd]
                if k in blsrcs: rsvdat[k][bl] = rsvdat[k].get(bl,[]) + [rsvd]
                blwgt[k][bl] = blwgt[k].get(bl,[]) + [_f * n.sqrt(u**2 + v**2)]
    ants = ants.keys(); ants.sort()
    # Simulate visibilities for each source for the data we collected
    for bl in msrdat:
        msrdat[bl] = n.array(msrdat[bl])
        msrval[bl] = n.array(msrval[bl])
        for k in simdat.keys():
            simdat[k][bl] = n.array(simdat[k][bl])
            if k in blsrcs: rsvdat[k][bl] = n.array(rsvdat[k][bl])
            blwgt[k][bl] = n.array(blwgt[k][bl])
            if n.all(simdat[k][bl] == 0):
                del(simdat[k])
                del(blwgt[k])
                if k in blsrcs: del(rsvdat[k])
    steps = [{'srcs':{}, 'score': n.sqrt(sum([n.sum(n.abs(msrdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))}]
    print 'Initial Score: %f' % (steps[0]['score'])

    # Initialize filters
    SHAPE = msrdat.values()[0].shape
    filter = gen_filter(SHAPE, opts.dw, opts.drw)
    filter_take = n.where(filter)
    window = gen_window(SHAPE)
    #def to_coeffs(d, noise=True):
    def to_coeffs(d, noise=False):
        c = ifft2(d*window)[filter_take]
        if noise:
            noise_amp = n.random.normal(scale=EXP_NOISE*d.size/n.sqrt(c.size), size=c.shape)
            #noise_amp = n.random.normal(scale=EXP_NOISE, size=c.shape)
            noise_phs = n.random.uniform(0,2*n.pi, size=c.shape)
            noise = noise_amp * n.exp(1j*noise_phs)
            #print n.sqrt(n.average(n.abs(noise)**2)),
            #print n.sqrt(n.average(n.abs(from_coeffs(noise))**2)),
            #print SHAPE, noise.shape
            c += noise
        return c

    def from_coeffs(c):
        d = n.zeros(SHAPE, dtype=n.complex)
        d[filter_take] = c
        return fft2(d) / d.size
        #return fft2(d)

    def compute_score(resdat, msrval):
        return n.sqrt(sum([n.sum(n.abs((resdat[bl]*window)**2)) for bl in msrval]) / (sum([n.sum(msrval[bl]*window) for bl in msrval])))

    # Create residuals
    resdat = {}
    for bl in msrdat: resdat[bl] = msrdat[bl].copy()

    dw, drw = 0, 0
    mask = coeff_mask(filter, dw, drw)
    for iter in xrange(opts.maxiter):
        print 'Iteration %d:' % iter
        # Derive source estimators
        step = {'srcs':{}, 'd_ant':{}}
        if iter == 0: # BEAMFORM
            for k in simdat:
                d, w = 0, 0
                for bl in simdat[k]:
                    if k in blsrcs: wgt = blwgt[k][bl] * rsvdat[k][bl]
                    else: wgt = blwgt[k][bl]
                    d += (resdat[bl] * n.conj(simdat[k][bl])).real * wgt
                    if k in blsrcs: wgt *= rsvdat[k][bl]
                    w += wgt
                d = n.where(w == 0, 0, n.sqrt(n.abs(d)/w))
                #c = to_coeffs(d, False)# * coeff_mask(filter, 0, 0)
                c = to_coeffs(d, False) * mask
                step['srcs'][k] = {}
                for i in ants: step['srcs'][k][i] = c.copy()
        else: # PER ANT
            # Sometimes change filter width
            if iter < BURN_IN:
                dw, drw = int(opts.dw*iter/BURN_IN), int(opts.drw*iter/BURN_IN)
                mask = coeff_mask(filter, dw, drw)
            elif iter == BURN_IN:
                dw,drw = opts.dw, opts.drw
                mask = 1
            print '    DW=%d DRW=%d' % (dw, drw)
                
            if iter < BURN_IN:
                for k in simdat:
                    d,w = {}, {}
                    for bl in simdat[k]:
                        i,j = a.miriad.bl2ij(bl)
                        if k in blsrcs: wgt = blwgt[k][bl] * rsvdat[k][bl]
                        else: wgt = blwgt[k][bl]
                        _d = resdat[bl] * n.conj(simdat[k][bl]) * wgt
                        if k in blsrcs: wgt *= rsvdat[k][bl]
                        _w = wgt
                        #d[i] = d.get(i,0) + .5*_d
                        d[i] = d.get(i,0) + _d
                        w[i] = w.get(i,0) + _w * n.conj(steps[-1]['ant'][k][j])
                        #d[j] = d.get(j,0) + n.conj(.5*_d)
                        d[j] = d.get(j,0) + n.conj(_d)
                        w[j] = w.get(j,0) + _w * n.conj(steps[-1]['ant'][k][i])
                    step['srcs'][k] = {}
                    step['d_ant'][k] = {}
                    for i in d:
                        # This divide can cause an explosion if w[i] is close to
                        # 0, which can happen early on when summing weights that
                        # depend on estimated antenna gain.  Adding noise helps
                        # This is also not exactly the formula we want: we should
                        # technically filter before dividing by w[i], but that
                        # results in a convolution in DDR domain...
                        #c = to_coeffs(n.where(w[i] == 0, 0, d[i]/w[i]))
                        step['d_ant'][k][i] = n.where(w[i] == 0, 0, d[i]/w[i])
                        c = to_coeffs(step['d_ant'][k][i])
                        #if iter < BURN_IN: c *= CLEAN_GAIN
                        if True: c *= CLEAN_GAIN
                        #step['srcs'][k][i] = (steps[-1]['srcs'][k][i] + c) #* _mask[k]
                        step['srcs'][k][i] = (steps[-1]['srcs'][k][i] + c) * mask
            # Estimate noise
            #if iter >= BURN_IN:
            if True:
                for k in simdat:
                    if not step['srcs'].has_key(k): step['srcs'][k] = {}
                    for i in d:
                        if not step['srcs'][k].has_key(i): step['srcs'][k][i] = steps[-1]['srcs'][k][i].copy()
                        #noise = 0
                        #for K in simdat:
                        #    _n,_w = 0.,0.
                        #    if k == K: continue
                        #    for j in d:
                        #        if i == j: continue
                        #        #_n += n.sum(n.abs(step['d_ant'][K][i]*n.conj(steps[-1]['ant'][K][j]) + steps[-1]['ant'][K][i]*n.conj(step['d_ant'][K][j]))**2)
                        #        #_n += n.sum((step['d_ant'][K][i]*n.conj(steps[-1]['ant'][K][j]) + steps[-1]['ant'][K][i]*n.conj(step['d_ant'][K][j])) * simdat[K][a.miriad.ij2bl(i,j)])
                        #        _n += n.sum(n.abs(step['d_ant'][K][i]*n.conj(steps[-1]['ant'][K][j]) + steps[-1]['ant'][K][i]*n.conj(step['d_ant'][K][j]))**2)
                        #        _w += steps[-1]['ant'][K][i].size
                        #        #print K, j, _n, _w
                        #    noise += _n/_w
                        ##noise = n.abs(noise)
                        #noise = n.sqrt(noise)
                        #print k, i, n.abs(step['srcs'][k][i][0]), n.abs(steps[-1]['srcs'][k][i][0]), noise,
                        #noise = steps[-1]['score'] / n.sqrt(float(len(d)))
                        #noise = steps[-1]['score']
                        #noise_amp = n.random.normal(scale=noise, size=step['srcs'][k][i].shape)
                        #print k, i, n.abs(step['srcs'][k][i])
                        noise = n.abs(step['srcs'][k][i])
                        noise_amp = n.random.normal(scale=FRAC_NOISE*noise.clip(NOISE_CLIP*noise.max(),n.Inf))
                        noise_phs = n.random.uniform(0,2*n.pi, size=step['srcs'][k][i].shape)
                        noise = noise_amp * n.exp(1j*noise_phs)
                        step['srcs'][k][i] += noise
                        #print n.abs(step['srcs'][k][i][0])
        # Compute residual for all sources
        _resdat = compute_residual(msrdat, msrval, simdat, rsvdat, step)
        step['score'] = compute_score(_resdat,msrval)

        # Always accept any downhill jump, advance step if bottoming out
        #if n.random.uniform() < n.exp((steps[-1]['score'] - step['score'])/EXP_NOISE):
        #_deg_of_freedom = sum([_mask[k].sum() for k in _mask])
        #deg_of_freedom = sum([mask[k].sum() for k in mask])
        #print deg_of_freedom, _deg_of_freedom, 1.01**(deg_of_freedom - _deg_of_freedom),
        print '    Probability:', n.exp((steps[-1]['score'] - step['score'])/EXP_NOISE)
        #if n.random.uniform() < n.exp((steps[-1]['score'] - step['score'])/EXP_NOISE) * 1.01**(deg_of_freedom - _deg_of_freedom):
        #if iter < BURN_IN or n.random.uniform() < n.exp((steps[-1]['score'] - step['score'])/EXP_NOISE):
        if n.random.uniform() < n.exp((steps[-1]['score'] - step['score'])/EXP_NOISE):
        #if True:
            if len(steps) > 1:
                del(steps[-1]['ant'])
                del(steps[-1]['xtalk'])
            steps.append(step)
            resdat = _resdat
            #dw, drw = _dw.copy(), _drw.copy()
            #mask = _mask.copy()
            h, bins = n.histogram([s['score'] for s in steps], bins=500, range=(0,1000))
            h = n.log10(h.astype(n.float) / h.max()).clip(-3,0)
            coeffs = {}
            h2 = {}
            #for k in step['srcs']:
            #for k in ['cyg','cas','21:44:17.81_28:12:24.7']:
            for k in ['cyg','cas']:
                coeffs[k] = n.array([s['srcs'][k][i] for s in steps[1:] for i in s['srcs'][k]])
                _c = coeffs[k].flatten()
                _c_abs = n.abs(_c)
                #_c *= (n.log10(_c_abs).clip(1,6)-1) / _c_abs
                _c *= (n.log10(_c_abs).clip(.5,5.5)-.5) / _c_abs
                h2[k],xedges,yedges = n.histogram2d(_c.real, _c.imag, bins=[400,400], range=[[-5,5],[-5,5]])
                h2[k] = n.log10(h2[k].astype(n.float) / h2[k].max()).clip(-3,0)
            if LINE1 is None:
                P.subplot(231)
                bins = .5*(bins[1:] + bins[:-1])
                LINE1, = P.plot(bins, h, 'k-')
                P.xlim(0,1000)
                P.ylim(-3,0)
                for i,k in enumerate(coeffs):
                    P.subplot(233)
                    LINE[k], = P.semilogy(n.average(n.abs(from_coeffs(n.average(coeffs[k], axis=0)))**2, axis=1))
                    P.ylim(1e-1,1e5)
                    P.subplot(2,3,4+i)
                    LINE2[k] = P.imshow(h2[k], origin='lower', aspect='auto', vmin=-3, vmax=0, extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]))
            else:
                LINE1.set_ydata(h)
                for k in coeffs:
                    LINE[k].set_ydata(n.average(n.abs(from_coeffs(n.average(coeffs[k], axis=0)))**2, axis=1))
                    LINE2[k].set_data(h2[k])
        else:
            print '    Rejecting proposal: Score %f (%5.2f%%)' % (step['score'], n.round(100*step['score']/steps[0]['score'],2))

        P.draw()
        print '    New Score: %f (%5.2f%%)' % (steps[-1]['score'], n.round(100*steps[-1]['score']/steps[0]['score'],2))
    print '    Final Score: %f (%5.2f%%)' % (steps[-1]['score'], n.round(100*steps[-1]['score']/steps[0]['score'],2))

    # Pare down data and save to file
    if TRIM:
        i0 = len(times) / 6
        times = times[i0:-i0]
    n.savez('%s__times.npz' % (arg), times=n.array(times))
    n.savez('%s__afreqs.npz' % (arg), freqs=afreqs)
    __xtalk = {}
    for k in step['srcs']:
        d = {}
        for i in step['srcs'][k]:
            if TRIM: d[str(i)] = step['srcs'][k][i][i0:-i0]
            else: d[str(i)] = step['srcs'][k][i]
        n.savez( '%s__step__%s.npz' % (arg,k), **d)

