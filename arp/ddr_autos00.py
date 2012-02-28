#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, math, glob
PLOT = True
if PLOT: import pylab as P; P.ion()

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('-r', '--drw', dest='drw', type=int, default=5,
    help='The number of delay-rate bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('--maxiter', dest='maxiter', type='int', default=100,
    help='Maximum number of iterations to allow.')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])

MAX_CASCADE_LEN = max(opts.dw,opts.drw) * 3
#MAX_CASCADE_LEN = 20
MIN_CASCADE_LEN = MAX_CASCADE_LEN / 3
CASCADE_PROBABILITY = .005
CLEAN_GAIN = 0.5
#NOISE_CLIP = 1e-2
NOISE_CLIP = 10**0.5
#EXP_NOISE = 5.
#FRAC_NOISE = n.sqrt(EXP_NOISE) * 1e-3
#FRAC_NOISE = 2e-3
#EXP_NOISE = FRAC_NOISE * 1e3 * .75
FRAC_NOISE = 1e-3
#EXP_NOISE = FRAC_NOISE * 1e3 * .025
EXP_NOISE = .005

def rms(d): return n.sqrt(n.average(n.abs(d)**2))
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
    x,y = n.indices(shape)
    x = n.cos(2*n.pi/shape[0]*x)
    y = n.cos(2*n.pi/shape[1]*y)
    w = (1 - x) * (1 - y)
    return w

try:
    import fftw3, cow
    print 'Using FFTW FFT'
    _fftw3_dat, _fftw3_fwd, _fftw3_rev = None, None, None
    def fft2(d):
        global _fftw3_fwd, _fftw3_dat
        if _fftw3_fwd is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_fwd = fftw3.Plan(_fftw3_dat, None, direction='forward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_fwd()
        return _fftw3_dat.copy()
    def ifft2(d):
        global _fftw3_rev, _fftw3_dat
        if _fftw3_rev is None:
            if _fftw3_dat is None: _fftw3_dat = n.zeros(d.shape, dtype=n.complex)
            _fftw3_rev = fftw3.Plan(_fftw3_dat, None, direction='backward', flags=['measure'])
        _fftw3_dat[:] = d
        _fftw3_rev()
        return _fftw3_dat.copy()
except(ImportError): 
    print 'Using numpy FFT'
    fft2, ifft2 = n.fft.fft2, n.fft.ifft2

def compute_residual(msrdat, msrval, step):
    # Compute residual for all sources
    _resdat = {}
    for bl in msrdat:
        i,j = a.miriad.bl2ij(bl)
        _resdat[bl] = msrdat[bl] - from_coeffs(step['coeffs'][bl]) * msrval[bl]
    return _resdat


LINE1, LINE2 = None, None

for filename in args:
    # Gather data
    print 'Processing', filename
    times = []
    # Collect data
    msrdat, msrval = {}, {}
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or t != times[-1]:
            times.append(t)
        bl = a.miriad.ij2bl(i,j)
        f = n.logical_not(f).astype(n.float)
        d = d * f #/ aa.passband(i,j)
        d -= n.abs(d.sum() / max(f.sum(),1)) * f
        msrdat[bl] = msrdat.get(bl,[]) + [d]
        msrval[bl] = msrval.get(bl,[]) + [f]
    # Simulate visibilities for each source for the data we collected
    for bl in msrdat:
        msrdat[bl] = n.array(msrdat[bl])
        msrval[bl] = n.array(msrval[bl])
    # Initialize filters
    SHAPE = msrdat.values()[0].shape
    print SHAPE
    filter = gen_filter(SHAPE, opts.dw, opts.drw)
    filter_take = n.where(filter)
    window = gen_window(SHAPE)
    resdat, ker, coeffs = {}, {}, {}
    for bl in msrdat:
        msrdat[bl] *= window
        msrval[bl] *= window
        resdat[bl] = msrdat[bl].copy()
        k = ifft2(msrval[bl])
        g = rms(msrval[bl])
        ker[bl] = (k,g)
        coeffs[bl] = n.zeros((filter_take[0].size,), dtype=n.complex)

    def to_coeffs(d, bl):
        _d = ifft2(d)
        k,g = ker[bl]
        _d, info = a.deconv.clean(_d, k, tol=1e-3)
        _d += info['res'] / g
        return _d[filter_take]
    def from_coeffs(c):
        d = n.zeros(SHAPE, dtype=n.complex)
        d[filter_take] = c
        #return fft2(d) / d.size
        return fft2(d)
    #def compute_score(resdat, msrval):
    #    return n.sqrt(sum([n.sum(n.abs((resdat[bl]*window)**2)) for bl in msrval]) / (sum([n.sum(msrval[bl]*window) for bl in msrval])))
    def compute_score(resdat, msrval):
        return n.sqrt(sum([n.sum(n.abs((resdat[bl])**2)) for bl in msrval]) / (sum([n.sum(msrval[bl]) for bl in msrval])))

    step = {'coeffs':coeffs, 'burn':True}
    step['score'] = compute_score(msrdat,msrval)
    steps = [step]
    # Main Peeling Loop
    for iter in xrange(opts.maxiter):
        print 'Iteration %d:' % iter
        
        # Derive source estimators
        step = {'coeffs':{}, 'burn':False}
        prev = steps[-1]
        if iter == 0 or n.random.uniform() < CASCADE_PROBABILITY:
            # Backtrack some fraction for each source
            for bl in prev['coeffs']:
                F = n.random.uniform(size=prev['coeffs'][bl].shape)
                step['coeffs'][bl] = F * prev['coeffs'][bl]
            resdat = compute_residual(msrdat, msrval, step)
            step['score'] = compute_score(resdat,msrval)
            print '    Seed Score: %f' % (step['score'])
            cascade_len = n.random.randint(MIN_CASCADE_LEN, MAX_CASCADE_LEN)
            for i in xrange(cascade_len):
                _step = {'coeffs':{}, 'burn':False}
                for bl in step['coeffs']:
                    _step['coeffs'][bl] = step['coeffs'][bl].copy()
                print '    CASCADE %d/%d:' % (i+1,cascade_len),
                dw = int(float(opts.dw+1)*i/cascade_len) % (opts.dw+1)
                drw = int(float(opts.drw+1)*i/cascade_len) % (opts.drw+1)
                #dw = opts.dw
                #drw = opts.drw
                print 'DW=%d DRW=%d' % (dw, drw)
                mask = coeff_mask(filter, dw, drw)
                for bl in resdat:
                    G = max(0,n.random.uniform(-CLEAN_GAIN/2,CLEAN_GAIN))
                    #G = 1
                    c = to_coeffs(resdat[bl],bl)
                    _step['coeffs'][bl] += G * c * mask
                # Add noise
                for bl in _step['coeffs']:
                    noise = n.abs(_step['coeffs'][bl])
                    noise_amp = n.random.normal(scale=FRAC_NOISE*noise.clip(NOISE_CLIP,n.Inf))
                    noise_phs = n.random.uniform(0,2*n.pi, size=_step['coeffs'][bl].shape)
                    noise = noise_amp * n.exp(1j*noise_phs)
                    _step['coeffs'][bl] += noise
                # Compute residual for all sources
                _resdat = compute_residual(msrdat, msrval, _step)
                _step['score'] = compute_score(_resdat,msrval)

                print '        Prob:', n.exp((step['score'] - _step['score'])/EXP_NOISE)
                if n.random.uniform() < n.exp((step['score'] - _step['score'])/EXP_NOISE):
                    step = _step
                    resdat = _resdat
                else:
                    print '        Rejecting: Score %f (%5.2f%%)' % (_step['score'], n.round(100*_step['score']/steps[0]['score'],2))
                print '        Score: %f (%5.2f%%)' % (step['score'], n.round(100*step['score']/steps[0]['score'],2))
        else:
            for bl in prev['coeffs']:
                step['coeffs'][bl] = prev['coeffs'][bl].copy()
                
        for bl in resdat:
            G = max(0,n.random.uniform(-CLEAN_GAIN/2,CLEAN_GAIN))
            #G = 1
            c = to_coeffs(resdat[bl],bl)
            step['coeffs'][bl] += G * c
            
        # Add noise
        for bl in step['coeffs']:
            noise = n.abs(step['coeffs'][bl])
            noise_amp = n.random.normal(scale=FRAC_NOISE*noise.clip(NOISE_CLIP,n.Inf))
            noise_phs = n.random.uniform(0,2*n.pi, size=step['coeffs'][bl].shape)
            noise = noise_amp * n.exp(1j*noise_phs)
            step['coeffs'][bl] += noise
        # Compute residual for all sources
        _resdat = compute_residual(msrdat, msrval, step)
        step['score'] = compute_score(_resdat,msrval)

        # Always accept any downhill jump, advance step if bottoming out
        print '    Prob:', n.exp((prev['score'] - step['score'])/EXP_NOISE)
        if n.random.uniform() < n.exp((prev['score'] - step['score'])/EXP_NOISE):
            # Clear off unnecessary data
            steps.append(step)
            resdat = _resdat
            if PLOT:
                _d = sum([resdat[bl] for bl in resdat]) / len(resdat)
                plot_d1 = n.log10(n.abs(_d/window.clip(1e-3,n.Inf))).clip(-10,n.Inf)
                #plot_d1 = n.log10(n.abs(resdat[bl])).clip(-10,n.Inf)
                _d = ifft2(_d)
                k,g = ker[bl]
                _d, info = a.deconv.clean(_d, k, tol=1e-3)
                _d += info['res'] / g
                plot_d2 = n.log10(n.abs(_d)).clip(-10,n.Inf)
                plot_d2 = a.img.recenter(plot_d2, n.array(plot_d2.shape)/2)
                if LINE1 is None:
                    P.subplot(1,2,1)
                    #mx = plot_d1.max()
                    mx = 0
                    LINE1 = P.imshow(plot_d1, vmax=mx, vmin=mx-2, aspect='auto', interpolation='nearest')
                    #P.colorbar(shrink=.5)
                    P.subplot(1,2,2)
                    #mx = plot_d2.max()
                    mx = -3
                    LINE2 = P.imshow(plot_d2, vmax=mx, vmin=mx-2, aspect='auto', interpolation='nearest')
                else:
                    LINE1.set_data(plot_d1)
                    LINE2.set_data(plot_d2)
        else:
            print '    Rejecting: Score %f (%5.2f%%)' % (step['score'], n.round(100*step['score']/steps[0]['score'],2))
            # Clear off unnecessary data
            step = steps[-1].copy()
            steps.append(step)

        # Draw plots
        if PLOT: P.draw()
        print '    Score: %f (%5.2f%%)' % (steps[-1]['score'], n.round(100*steps[-1]['score']/steps[0]['score'],2))
    print 'Final Score: %f (%5.2f%%)' % (steps[-1]['score'], n.round(100*steps[-1]['score']/steps[0]['score'],2))

    # Pare down data and save to file
    good_steps = [s for s in steps if not s['burn']]
    scores = n.array([s['score'] for s in good_steps])
    n.savez('%s__autoinfo.npz' % (filename), times=n.array(times), chans=SHAPE[1], scores=scores)
    d = {}
    for bl in good_steps[-1]['coeffs']:
        d[str(bl)] = n.array([s['coeffs'][bl] for s in good_steps])
    n.savez( '%s__auto.npz' % (filename), **d)

#import time
#while True:
#    P.draw()
#    time.sleep(.01)
