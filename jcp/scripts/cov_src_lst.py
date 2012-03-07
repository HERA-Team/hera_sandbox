#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, math, glob, os
PLOT = False
if PLOT: import pylab as P; P.ion()

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, chan=True, src=True,
    dec=True, ant=True, pol=True)
o.add_option('-l', '--lstrng', dest='lstrng', default='0_6.2832',
    help='Range of LSTs to bin.')
o.add_option('-b', '--blsrcs', dest='blsrcs', default='',
    help='Sources (from any tier) that need per-baseline solutions')
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-r', '--drw', dest='drw', type=int, default=5,
    help='The number of delay-rate bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('--maxiter', dest='maxiter', type='int', default=100,
    help='Maximum number of iterations to allow.')
opts, args = o.parse_args(sys.argv[1:])

# Select only some of these LST bins to analyze
lstrng = map(float, opts.lstrng.split('_'))
    
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
NTIERS = len(srctier)
cat = a.fit.SrcCatalog()
for c in srctier: cat.add_srcs(c.values())
blsrcs = opts.blsrcs.split(',')

MAX_CASCADE_LEN = 200
MIN_CASCADE_LEN = 20
CASCADE_PROBABILITY = .002
#CASCADE_PROBABILITY = .0
CLEAN_GAIN = 0.5
#NOISE_CLIP = 1e-2
NOISE_CLIP = 10**0.5
#EXP_NOISE = 5.
#FRAC_NOISE = n.sqrt(EXP_NOISE) * 1e-3
FRAC_NOISE = 2e-3
EXP_NOISE = FRAC_NOISE * 1e3 * .25
#FRAC_NOISE = 1e-2
#EXP_NOISE = FRAC_NOISE * 1e3 * .025

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
    if not step.has_key('ant'): step['ant'] = {}
    for k in step['srcs']:
        if not step['ant'].has_key(k): step['ant'][k] = {}
        for i in step['srcs'][k]:
            if step['ant'][k].has_key(i): continue
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

LINE1, LINE, LINE2 = None, {}, {}
def update_plot(steps):
    global LINE, LINE1, LINE2
    step = steps[-1]
    good_steps = [s for s in steps if not s['burn']]
    if len(good_steps) < 10: return
    h, bins = n.histogram([s['score'] for s in good_steps], bins=500, range=(0,1000))
    h = n.log10(h.astype(n.float) / h.max()).clip(-3,0)
    coeffs = {}
    h2 = {}
    plot_srcs = [k for k in step['srcs'] if not k in ['cyg','cas']]
    #for k in ['cyg','cas'] + plot_srcs[:2]:
    for k in ['vir','cyg'] + plot_srcs[:2]:
        coeffs[k] = n.array([s['srcs'][k][i] for s in good_steps for i in s['srcs'][k]])
        _c = coeffs[k].flatten()
        _c_abs = n.abs(_c)
        _c *= (n.log10(_c_abs).clip(0.5,6)) / _c_abs
        h2[k],xedges,yedges = n.histogram2d(_c.real, _c.imag, bins=[500,500], range=[[-6,6],[-6,6]])
        h2[k] = n.log10(h2[k].astype(n.float) / h2[k].max()).clip(-3,0)
    if LINE1 is None:
        P.subplot(231)
        bins = .5*(bins[1:] + bins[:-1])
        LINE1, = P.plot(bins, h, 'k-')
        P.xlim(0,1000)
        P.ylim(-3,0)
        for pltnum,k in enumerate(coeffs):
            _c = n.average(n.average([n.abs(from_coeffs(step['srcs'][k][i]))**2 for i in step['srcs'][k]], axis=0), axis=1)
            P.subplot(232)
            LINE[k], = P.semilogy(_c)
            P.ylim(1e-1,1e5)
            P.subplot(2,3,3+pltnum)
            LINE2[k] = P.imshow(h2[k], origin='lower', aspect='auto', 
                vmin=-3, vmax=0, extent=(yedges[0],yedges[-1],xedges[0],xedges[-1]))
            P.title(k)
    else:
        LINE1.set_ydata(h)
        for k in coeffs:
            _c = n.average(n.average([n.abs(from_coeffs(step['srcs'][k][i]))**2 for i in step['srcs'][k]], axis=0), axis=1)
            LINE[k].set_ydata(_c)
            LINE2[k].set_data(h2[k])

file_glob = '.'.join(map(globize, args[0].split('.')))
filelist = glob.glob(file_glob)
filelist.sort()


gargs = []
for arg in args:
    uv = a.miriad.UV(arg)
    if lstrng[0] < lstrng[1]:
        uv.select('lst', lstrng[0], lstrng[1])
    else:
        uv.select('lst', lstrng[1], lstrng[0], include=False)
    try:
        p, d = uv.read()
        gargs.append(arg)
    except:
        continue

if len(gargs) == 0: print 'No files have that lst range.'

for arg in gargs:
    # Gather data
    print 'Processing', arg
    if os.path.exists('%s__info.npz' % (arg)):
        print '%s exists. Skipping...' % ('%s__info.npz' % (arg))
        continue
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
        if lstrng[0] < lstrng[1]:
            uv.select('lst', lstrng[0], lstrng[1])
        else:
            uv.select('lst', lstrng[1], lstrng[0], include=False)
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
                except(a.phs.PointingError):
                    _f = n.zeros_like(f)
                    simd = n.zeros_like(d)
                    if k in blsrcs: rsvd = simd
                simdat[k][bl] = simdat[k].get(bl,[]) + [simd]
                if k in blsrcs: rsvdat[k][bl] = rsvdat[k].get(bl,[]) + [rsvd]
                blwgt[k][bl] = blwgt[k].get(bl,[]) + [_f]
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
    # Initialize filters
    SHAPE = msrdat.values()[0].shape
    filter = gen_filter(SHAPE, opts.dw, opts.drw)
    filter_take = n.where(filter)
    window = gen_window(SHAPE)
    def to_coeffs(d):
        c = ifft2(d*window)[filter_take]
        return c
    def from_coeffs(c):
        d = n.zeros(SHAPE, dtype=n.complex)
        d[filter_take] = c
        return fft2(d) / d.size
        #return fft2(d)
    def compute_score(resdat, msrval):
        return n.sqrt(sum([n.sum(n.abs((resdat[bl]*window)**2)) for bl in msrval]) / (sum([n.sum(msrval[bl]*window) for bl in msrval])))

    steps = []

    # Compute initial score
    resdat = {}
    for bl in msrdat: resdat[bl] = msrdat[bl].copy()
    step = {'srcs':{}, 'burn':True}
    step['score'] = compute_score(resdat, msrval)
    steps.append(step)
    print 'Initial Score: %f' % (steps[-1]['score'])

    # Seed starting position
    step = {'srcs':{}, 'burn':True}
    mask = coeff_mask(filter, 0, 0)
    for k in simdat:
        d, w = 0, 0
        for bl in simdat[k]:
            if k in blsrcs: wgt = blwgt[k][bl] * rsvdat[k][bl]
            else: wgt = blwgt[k][bl]
            d += (resdat[bl] * n.conj(simdat[k][bl])).real * wgt
            if k in blsrcs: wgt *= rsvdat[k][bl]
            w += wgt
        d = n.where(w == 0, 0, n.sqrt(n.abs(d)/w))
        c = to_coeffs(d) * mask
        step['srcs'][k] = {}
        for i in ants: step['srcs'][k][i] = c.copy()
    resdat = compute_residual(msrdat, msrval, simdat, rsvdat, step)
    step['score'] = compute_score(resdat,msrval)
    steps.append(step)
    print 'Seed Score: %f' % (steps[-1]['score'])
    cascade_len = n.random.randint(MIN_CASCADE_LEN, MAX_CASCADE_LEN)

    # Main Peeling Loop
    for iter in xrange(opts.maxiter):
        print 'Iteration %d:' % iter
        # Derive source estimators
        step = {'srcs':{}, 'burn':False}
        prev = steps[-1]
        if iter == 0 or n.random.uniform() < CASCADE_PROBABILITY:
            # Backtrack some fraction for each source
            for k in prev['srcs']:
                step['srcs'][k] = {}
                for i in prev['srcs'][k]:
                    F = n.random.uniform(size=prev['srcs'][k][i].shape)
                    step['srcs'][k][i] = F * prev['srcs'][k][i]
            resdat = compute_residual(msrdat, msrval, simdat, rsvdat, step)
            step['score'] = compute_score(resdat,msrval)
            print '    Seed Score: %f' % (step['score'])
            cascade_len = n.random.randint(MIN_CASCADE_LEN, MAX_CASCADE_LEN)
            for i in xrange(cascade_len):
                _step = {'srcs':{}, 'burn':False}
                for k in step['srcs']:
                    _step['srcs'][k] = {}
                    for ai in step['srcs'][k]:
                        _step['srcs'][k][ai] = step['srcs'][k][ai].copy()
                print '    CASCADE %d/%d:' % (i+1,cascade_len),
                dw = int(float(opts.dw+1)*i*NTIERS/cascade_len) % (opts.dw+1)
                drw = int(float(opts.drw+1)*i*NTIERS/cascade_len) % (opts.drw+1)
                tier = NTIERS * i / cascade_len
                print 'DW=%d DRW=%d TIER=%d' % (dw, drw, tier)
                mask = coeff_mask(filter, dw, drw)
                compute_ant(step)
                for k in simdat:
                    G = max(0,n.random.uniform(-CLEAN_GAIN/2,CLEAN_GAIN))
                    if G == 0 or not k in srctier[tier]: continue
                    #print '        GAIN: %s %f' % (k, G)
                    d,w = {}, {}
                    phs,amp = {}, {}
                    for bl in simdat[k]:
                        i,j = a.miriad.bl2ij(bl)
                        if k in blsrcs: wgt = blwgt[k][bl] * rsvdat[k][bl]
                        else: wgt = blwgt[k][bl]
                        _d = resdat[bl] * n.conj(simdat[k][bl]) * wgt
                        if k in blsrcs: wgt *= rsvdat[k][bl]
                        _w = wgt
                        if not phs.has_key(i):
                            amp[i] = n.abs(step['ant'][k][i])
                            phs[i] = step['ant'][k][i] / amp[i]
                        if not phs.has_key(j):
                            amp[j] = n.abs(step['ant'][k][j])
                            phs[j] = step['ant'][k][j] / amp[j]
                        #d[i] = d.get(i,0) + _d
                        d[i] = d.get(i,0) + _d * phs[j]
                        #w[i] = w.get(i,0) + _w * n.conj(step['ant'][k][j])
                        w[i] = w.get(i,0) + _w * amp[j]
                        #d[j] = d.get(j,0) + n.conj(_d)
                        d[j] = d.get(j,0) + n.conj(_d) * phs[i]
                        #w[j] = w.get(j,0) + _w * n.conj(step['ant'][k][i])
                        w[j] = w.get(j,0) + _w * amp[i]
                    for i in d:
                        # This divide can cause an explosion if w[i] is close to
                        # 0, which can happen early on when summing weights that
                        # depend on estimated antenna gain.  Adding noise helps
                        # This is also not exactly the formula we want: should
                        # technically filter before dividing by w[i], but that
                        # results in a convolution in DDR domain...
                        c = to_coeffs(n.where(w[i] == 0, 0, d[i]/w[i]))
                        _step['srcs'][k][i] += G * c * mask
                # Add noise
                for k in _step['srcs']:
                    for i in _step['srcs'][k]:
                        noise = n.abs(_step['srcs'][k][i])
                        noise_amp = n.random.normal(scale=FRAC_NOISE*noise.clip(NOISE_CLIP,n.Inf))
                        noise_phs = n.random.uniform(0,2*n.pi, size=_step['srcs'][k][i].shape)
                        noise = noise_amp * n.exp(1j*noise_phs)
                        _step['srcs'][k][i] += noise
                # Compute residual for all sources
                _resdat = compute_residual(msrdat, msrval, simdat, rsvdat, _step)
                _step['score'] = compute_score(_resdat,msrval)

                print '        Prob:', n.exp((step['score'] - _step['score'])/EXP_NOISE)
                if n.random.uniform() < n.exp((step['score'] - _step['score'])/EXP_NOISE):
                    step = _step
                    resdat = _resdat
                else:
                    print '        Rejecting: Score %f (%5.2f%%)' % (_step['score'], n.round(100*_step['score']/steps[0]['score'],2))
                print '        Score: %f (%5.2f%%)' % (step['score'], n.round(100*step['score']/steps[0]['score'],2))
        else:
            for k in prev['srcs']:
                step['srcs'][k] = {}
                for i in prev['srcs'][k]:
                    step['srcs'][k][i] = prev['srcs'][k][i].copy()
                
        #print '    DW=%d DRW=%d TIER=%d' % (dw, drw, tier)
            
        # Add noise
        for k in step['srcs']:
            for i in step['srcs'][k]:
                noise = n.abs(step['srcs'][k][i])
                noise_amp = n.random.normal(scale=FRAC_NOISE*noise.clip(NOISE_CLIP,n.Inf))
                noise_phs = n.random.uniform(0,2*n.pi, size=step['srcs'][k][i].shape)
                noise = noise_amp * n.exp(1j*noise_phs)
                step['srcs'][k][i] += noise
        # Compute residual for all sources
        _resdat = compute_residual(msrdat, msrval, simdat, rsvdat, step)
        step['score'] = compute_score(_resdat,msrval)

        # Always accept any downhill jump, advance step if bottoming out
        print '    Prob:', n.exp((prev['score'] - step['score'])/EXP_NOISE)
        if n.random.uniform() < n.exp((prev['score'] - step['score'])/EXP_NOISE):
            # Clear off unnecessary data
            steps.append(step)
            resdat = _resdat
        else:
            print '    Rejecting: Score %f (%5.2f%%)' % (step['score'], n.round(100*step['score']/steps[0]['score'],2))
            # Clear off unnecessary data
            step = steps[-1].copy()
            steps.append(step)
        try:
            del(steps[-2]['ant'])
            del(steps[-2]['xtalk'])
        except(KeyError): pass

        # Draw plots
        if PLOT:
            if iter % 10 == 0: update_plot(steps)
            P.draw()
        print '    Score: %f (%5.2f%%)' % (steps[-1]['score'], n.round(100*steps[-1]['score']/steps[0]['score'],2))
    print 'Final Score: %f (%5.2f%%)' % (steps[-1]['score'], n.round(100*steps[-1]['score']/steps[0]['score'],2))

    # Pare down data and save to file
    good_steps = [s for s in steps if not s['burn']]
    scores = n.array([s['score'] for s in good_steps])
    n.savez('%s__info.npz' % (arg), times=n.array(times), freqs=afreqs, scores=scores)
    for k in good_steps[-1]['srcs']:
        d = {}
        for i in good_steps[-1]['srcs'][k]:
            d[str(i)] = n.array([s['srcs'][k][i] for s in good_steps])
        n.savez( '%s__%s.npz' % (arg,k), **d)

#import time
#while True:
#    P.draw()
#    time.sleep(.01)
