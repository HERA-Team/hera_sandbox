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

MAX_CASCADE_LEN = max(opts.dw,opts.drw) * 1
MAX_CASCADE_LEN = 10
MIN_CASCADE_LEN = MAX_CASCADE_LEN / 3
CASCADE_PROBABILITY = .005
CLEAN_GAIN = 0.5
NOISE_CLIP = 1e-10
#NOISE_CLIP = 10**0.5
#EXP_NOISE = 5.
#FRAC_NOISE = n.sqrt(EXP_NOISE) * 1e-3
#FRAC_NOISE = 2e-3
#EXP_NOISE = FRAC_NOISE * 1e3 * .75
#FRAC_NOISE = 1e-3
FRAC_NOISE = 1e-5
#EXP_NOISE = FRAC_NOISE * 1e3 * .025
EXP_NOISE = .0001

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
    #w = (1 - x) * (1 - y)
    w = (1 - y)
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
    return msrdat - from_coeffs(step['coeffs']) * msrval


LINE1, LINE2 = None, None
LINES = {}

dat, val = 0, 0
spec_sum, spec_wgt = 0, 0
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
        spec_sum += d
        spec_wgt += f
        #print d.compress(f).size, n.median(d.compress(f)), n.abs(d.sum() / max(f.sum(),1))
        #d -= n.abs(d.sum() / max(f.sum(),1)) * f
        if not n.all(f == 0): d /= n.abs(n.median(d[550:650].compress(f[550:650]))).clip(1e-3,n.Inf)
        msrdat[bl] = msrdat.get(bl,[]) + [d]
        msrval[bl] = msrval.get(bl,[]) + [f]
    # Simulate visibilities for each source for the data we collected
    for bl in msrdat:
        dat += n.array(msrdat[bl])
        val += n.array(msrval[bl])
spec_avg = spec_sum / spec_wgt.clip(1,n.Inf)
msrdat, msrval = dat, val
# Initialize filters
SHAPE = msrdat.shape
print SHAPE
filter = gen_filter(SHAPE, opts.dw, opts.drw)
filter_take = n.where(filter)
#window = gen_window(SHAPE)
resdat, ker = {}, {}
#msrdat *= window
#msrval *= window
resdat = msrdat.copy()
k = ifft2(msrval)
g = rms(msrval)
coeffs = n.zeros((filter_take[0].size,), dtype=n.complex)

def to_coeffs(d):
    _d = ifft2(d)
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
    return n.sqrt(n.sum(n.abs(resdat)**2) / n.sum(msrval**2))

step = {'coeffs':coeffs, 'burn':True}
step['score'] = compute_score(msrdat,msrval)
steps = [step]
# Main Peeling Loop
for iter in xrange(opts.maxiter):
    print 'Iteration %d:' % iter
    
    # Derive source estimators
    step = {'burn':False}
    prev = steps[-1]
    if iter == 0 or n.random.uniform() < CASCADE_PROBABILITY:
        # Backtrack some fraction for each source
        F = n.random.uniform(size=prev['coeffs'].shape)
        step['coeffs'] = F * prev['coeffs']
        resdat = compute_residual(msrdat, msrval, step)
        step['score'] = compute_score(resdat,msrval)
        print '    Seed Score: %f' % (step['score'])
        cascade_len = n.random.randint(MIN_CASCADE_LEN, MAX_CASCADE_LEN)
        for i in xrange(cascade_len):
            _step = {'burn':False}
            _step['coeffs'] = step['coeffs'].copy()
            print '    CASCADE %d/%d:' % (i+1,cascade_len),
            dw = int(float(opts.dw+1)*i/cascade_len) % (opts.dw+1)
            drw = int(float(opts.drw+1)*i/cascade_len) % (opts.drw+1)
            print 'DW=%d DRW=%d' % (dw, drw)
            mask = coeff_mask(filter, dw, drw)
            G = max(0,n.random.uniform(-CLEAN_GAIN/2,CLEAN_GAIN))
            c = to_coeffs(resdat)
            _step['coeffs'] += G * c * mask
            # Add noise
            noise = n.abs(_step['coeffs'])
            noise_amp = n.random.normal(scale=FRAC_NOISE*noise.clip(NOISE_CLIP,n.Inf))
            noise_phs = n.random.uniform(0,2*n.pi, size=_step['coeffs'].shape)
            noise = noise_amp * n.exp(1j*noise_phs)
            _step['coeffs'] += noise
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
        step['coeffs'] = prev['coeffs'].copy()
            
    G = max(0,n.random.uniform(-CLEAN_GAIN/2,CLEAN_GAIN))
    #G = 1
    c = to_coeffs(resdat)
    step['coeffs'] += G * c
        
    # Add noise
    noise = n.abs(step['coeffs'])
    noise_amp = n.random.normal(scale=FRAC_NOISE*noise.clip(NOISE_CLIP,n.Inf))
    noise_phs = n.random.uniform(0,2*n.pi, size=step['coeffs'].shape)
    noise = noise_amp * n.exp(1j*noise_phs)
    step['coeffs'] += noise
    # Compute residual for all sources
    sig = n.sqrt(n.sum(n.abs(_resdat*msrval)**2)/n.sum(n.abs(msrval)**2))
    goodval = n.where(n.abs(_resdat) < sig, 1, 0)
    _resdat = compute_residual(msrdat, msrval, step)
    step['score'] = compute_score(_resdat*goodval,msrval*goodval)

    # Always accept any downhill jump, advance step if bottoming out
    print '    Prob:', n.exp((prev['score'] - step['score'])/EXP_NOISE)
    if n.random.uniform() < n.exp((prev['score'] - step['score'])/EXP_NOISE):
        # Clear off unnecessary data
        steps.append(step)
        resdat = _resdat
        if PLOT:
            #plot_d1 = n.log10(n.abs(resdat/window.clip(1e-3,n.Inf))).clip(-10,n.Inf)
            #plot_d1 = n.log10(n.abs(resdat*goodval)).clip(-10,n.Inf)
            #plot_d1 = (resdat*goodval).real[:,550:650].flatten()
            plot_d1 = (resdat*goodval).real[:,592:594].flatten()
            print n.average(plot_d1), n.std(plot_d1), plot_d1.size
            h,bins = n.histogram(plot_d1.compress(n.abs(plot_d1) > 0))
            bins = (bins[1:] + bins[:-1])/2
            _d = ifft2(resdat*goodval)
            _d, info = a.deconv.clean(_d, k, tol=1e-3)
            _d += info['res'] #/ g
            __d = a.img.recenter(_d, n.array(_d.shape)/2)
            plot_d2 = n.log10(n.abs(__d)).clip(-10,n.Inf)
            plot_d3a = n.abs(__d[SHAPE[0]/2])
            plot_d3b = n.abs(__d[SHAPE[0]/2+1])
            plot_d3c = n.abs(__d[SHAPE[0]/2-1])
            #_d *= filter
            #_d[0] = 0
            #plot_d4a = n.abs(n.average(fft2(_d), axis=0))
            #plot_d4a = n.abs(n.sum(resdat*msrval, axis=0) / n.sum(msrval, axis=0))
            #plot_d4a = n.abs(n.sum((resdat-fft2(_d))*msrval*goodval, axis=0) / n.sum(msrval*goodval, axis=0))
            #plot_d4a.shape = (plot_d4a.size/8, 8)
            #plot_d4a = n.average(plot_d4a, axis=1)
            plot_d4b = n.abs(n.sum(resdat*msrval*goodval, axis=0) / n.sum(msrval*goodval, axis=0))
            #plot_d4b = n.abs(n.average(from_coeffs(step['coeffs']), axis=0))
            #plot_d4b.shape = (plot_d4b.size/8, 8)
            #plot_d4b = n.average(plot_d4b, axis=1)
            plot_d4c = n.abs(spec_avg)
            #plot_d4c.shape = (plot_d4c.size/8, 8)
            #plot_d4c = plot_d4c[:,0]
            if LINE1 is None:
                P.subplot(2,2,1)
                #mx = plot_d1.max()
                mx = -2
                #LINE1 = P.imshow(plot_d1, vmax=mx, vmin=mx-2, aspect='auto', interpolation='nearest')
                #LINE1 = P.imshow(plot_d1, vmax=1e-3, vmin=-1e-3, aspect='auto', interpolation='nearest')
                #P.colorbar(shrink=.5)
                LINE1, = P.plot(bins, h)
                P.subplot(2,2,2)
                #mx = plot_d2.max()
                mx = -4
                LINE2 = P.imshow(plot_d2, vmax=mx, vmin=mx-2, aspect='auto', interpolation='nearest')
                P.colorbar(shrink=.5)
                P.subplot(2,2,3)
                LINES[0], = P.semilogy(n.abs(plot_d3a))
                LINES[1], = P.semilogy(n.abs(plot_d3b))
                LINES[2], = P.semilogy(n.abs(plot_d3c))
                P.subplot(2,2,4)
                #LINES[10], = P.semilogy(n.abs(plot_d4a))
                LINES[11], = P.semilogy(n.abs(plot_d4b))
                LINES[12], = P.semilogy(n.abs(plot_d4c))
                P.ylim(1e-5,1e1)
            else:
                #LINE1.set_data(plot_d1)
                P.subplot(2,2,1); P.xlim(bins[0], bins[-1])
                LINE1.set_xdata(bins)
                LINE1.set_ydata(h)
                LINE2.set_data(plot_d2)
                LINES[0].set_ydata(plot_d3a)
                LINES[1].set_ydata(plot_d3b)
                LINES[2].set_ydata(plot_d3c)
                #LINES[10].set_ydata(plot_d4a)
                LINES[11].set_ydata(plot_d4b)
                LINES[12].set_ydata(plot_d4c)
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
d = n.array([s['coeffs'] for s in good_steps])
n.savez( '%s__auto.npz' % (filename), coeffs=d)

#import time
#while True:
#    P.draw()
#    time.sleep(.01)
