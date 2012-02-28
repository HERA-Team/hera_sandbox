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

CASCADE_LEN = 10
CLEAN_GAIN = 0.9
NOISE_CLIP = 1e-10
FRAC_NOISE = 1e-5
EXP_NOISE = .0001

def Tsys(chan, sdf=9.765e-05, sfreq=.1, t_sync=240, index=-2.5, t_rx=150.):
    return t_sync * ((chan * sdf + sfreq)/ .150)**index + t_rx

def Teor(chan, transition_chan=600, t_eor=.01):
    eor = t_eor * n.ones_like(chan)
    eor[transition_chan:] = 0
    return eor

bp_poly = [3708527239878.0649, -4532011120582.3359, 2417448183275.3481, -735123189051.80493, 139377660668.0575, -16870609639.249588, 1273053690.1204848, -54751484.793454409, 1027464.8904983159]

def rms(d,wgt=None):
    if wgt == None: return n.sqrt(n.average(n.abs(d)**2))
    else: return n.sqrt(n.sum(n.abs(d)**2) / n.sum(n.abs(wgt)**2))

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

def coeff_mask(filter, shape, red_dw, red_drw):
    filter_take = filter
    red_filter = gen_filter(shape, red_dw, red_drw)
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

fft2, ifft2 = n.fft.fft2, n.fft.ifft2

LINE1, LINE2 = None, None
LINES = {}

for filename in args:
    # Gather data
    print 'Reading', filename
    msrdat, msrval = {}, {}
    spec_sum, spec_wgt = 0, 0
    filters, kernels, gains = {}, {}, {}
    coeffs = {}
    times = []
    # Collect data
    uv = a.miriad.UV(filename)
    chan = n.arange(uv['nchan'])
    tsys = Tsys(chan, uv['sdf'], uv['sfreq'])
    teor = Teor(chan)
    bp = n.polyval(bp_poly, chan*uv['sdf']+uv['sfreq'])
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or t != times[-1]: times.append(t)
        bl = a.miriad.ij2bl(i,j)
        v = n.logical_not(f).astype(n.float)
        #print n.abs(n.median(d[550:650].compress(v[550:650]))).clip(1e-3,n.Inf)
        #if not n.all(v[550:650] == 0): d /= n.abs(n.median(d[550:650].compress(v[550:650]))).clip(1e-3,n.Inf)
        #d = bp * (tsys + teor)
        #d = tsys.copy()
        if not n.all(v[550:650] == 0): d /= n.abs(n.median(d[550:650].compress(v[550:650]))).clip(1e-3,n.Inf)
        else:
            d *= 0
            v *= 0
        d,v = d[60:-100], v[60:-100]
        d *= v
        msrdat[bl] = msrdat.get(bl,[]) + [d]
        msrval[bl] = msrval.get(bl,[]) + [v]
    chan = chan[60:-100]
    tsys = tsys[60:-100]
    #N = d.size - 160
    N = d.size + 100
    window = 0.35875 - 0.48829 * n.cos(2*n.pi*n.arange(N)/N) + 0.14128 * n.cos(4*n.pi*n.arange(N)/N) - 0.01168 * n.cos(6*n.pi*n.arange(N)/N)
    #window = n.concatenate([n.zeros(60), window, n.zeros(100)]) * 10
    window = window[50:-50] * 10
    #window *= 10
    window.shape = (1,window.size)
    for bl in msrdat:
        msrdat[bl] = n.array(msrdat[bl]) * window
        msrval[bl] = n.array(msrval[bl]) * n.where(window == 1, 0, 1)
        #msrval[bl][:,:70] = 0
        #msrval[bl][:,915:] = 0
    shape = msrval[bl].shape
    f = n.where(gen_filter(shape, opts.dw, opts.drw))
    #k = ifft2(msrval[bl]*window)
    #g = rms(msrval[bl]*window)
    k = ifft2(msrval[bl])
    g = rms(msrval[bl])
    for bl in msrdat: coeffs[bl] = n.zeros((f[0].size,), dtype=n.complex)

    def to_coeffs(d, shared=False, clean=1e-2):
        c = {}
        if shared:
            bl = d.keys()[0]
            _sum = d[bl]
            _wgt = msrval[bl]
            _d = ifft2(_sum/_wgt.clip(1,n.Inf))
            if clean > 0:
                _d, info = a.deconv.clean(_d, k, tol=clean)
                _d += info['res'] / g
            _c = _d[f]
            for bl in d: c[bl] = _c.copy()
        else:
            for bl in d:
                _d = ifft2(d[bl])
                if clean > 0:
                    _d, info = a.deconv.clean(_d, k, tol=clean)
                    _d += info['res'] / g
                c[bl] = _d[f]
        return c

    def from_coeffs(c):
        d = {}
        for bl in c:
            _d = n.zeros(msrval[bl].shape, dtype=n.complex)
            _d[f] = c[bl]
            d[bl] = fft2(_d) * msrval[bl]
        return d

    def compute_score(resdat, msrval):
        scr,wgt = 0,0
        for bl in msrval:
            scr += n.sum(n.abs(resdat[bl])**2)
            wgt += n.sum(n.abs(msrval[bl])**2)
        return n.sqrt(scr/wgt)

    def compute_residual(msrdat, msrval, step):
        # Compute residual for all sources
        d = from_coeffs(step['coeffs'])
        for bl in d: d[bl] = msrdat[bl] - d[bl]
        return d

    def deepcopy(d):
        _d = {}
        for bl in d: _d[bl] = d[bl].copy()
        return _d

    step = {'coeffs':coeffs, 'burn':True}
    step['score'] = compute_score(msrdat,msrval)
    steps = [step]

    # Main Peeling Loop
    for iter in xrange(opts.maxiter):
        print 'Iteration %d:' % iter
        
        # Derive source estimators
        step = {'burn':False}
        prev = steps[-1]
        if iter == 0:
            step['coeffs'] = deepcopy(prev['coeffs'])
            resdat = compute_residual(msrdat, msrval, step)
            step['score'] = compute_score(resdat,msrval)
            print '    Seed Score: %f' % (step['score'])
            for i in xrange(CASCADE_LEN):
                _step = {'burn':False}
                _step['coeffs'] = deepcopy(step['coeffs'])
                print '    CASCADE %d/%d:' % (i+1,CASCADE_LEN),
                dw = int(float(opts.dw+1)*i/CASCADE_LEN) % (opts.dw+1)
                drw = int(float(opts.drw+1)*i/CASCADE_LEN) % (opts.drw+1)
                print 'DW=%d DRW=%d' % (dw, drw)
                #c = to_coeffs(resdat)
                c = to_coeffs(resdat, shared=(i<CASCADE_LEN/2), clean=1e-3)
                G = CLEAN_GAIN
                mask = coeff_mask(f, msrval.values()[0].shape, dw, drw)
                for bl in _step['coeffs']: _step['coeffs'][bl] += G * c[bl] * mask
                # Compute residual for all sources
                _resdat = compute_residual(msrdat, msrval, _step)
                # Compute sig with old residual to avoid weird self-selections
                goodval = {}
                for bl in msrval:
                    #sig = n.sqrt(n.sum(n.abs(resdat[bl]*msrval[bl])**2)/n.sum(n.abs(msrval[bl])**2))
                    #goodval[bl] = n.where(n.abs(resdat[bl]) < sig, 1, 0) * msrval[bl]
                    goodval[bl] = msrval[bl]
                __resdat = {}
                for bl in _resdat: __resdat[bl] = _resdat[bl] * goodval[bl]
                _step['score'] = compute_score(__resdat,goodval)

                #_step['score'] = compute_score(_resdat,msrval)

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
        try: c = to_coeffs(__resdat, clean=0)
        except(NameError): c = to_coeffs(resdat, clean=0)
        for bl in step['coeffs']:
            step['coeffs'][bl] += G * c[bl]
            # Add noise
            noise = n.abs(step['coeffs'][bl])
            noise_amp = n.random.normal(scale=FRAC_NOISE*noise.clip(NOISE_CLIP,n.Inf))
            noise_phs = n.random.uniform(0,2*n.pi, size=step['coeffs'][bl].shape)
            noise = noise_amp * n.exp(1j*noise_phs)
            step['coeffs'][bl] += noise
        # Compute residual for all sources
        _resdat = compute_residual(msrdat, msrval, step)
        # Compute sig with old residual to avoid weird self-selections
        goodval = {}
        for bl in msrval:
            #sig = n.sqrt(n.sum(n.abs(resdat[bl]*msrval[bl])**2)/n.sum(n.abs(msrval[bl])**2))
            #goodval[bl] = n.where(n.abs(resdat[bl]) < sig, 1, 0) * msrval[bl]
            goodval[bl] = msrval[bl]
        __resdat = {}
        for bl in _resdat: __resdat[bl] = _resdat[bl] * goodval[bl]
        step['score'] = compute_score(__resdat,goodval)

        # Always accept any downhill jump, advance step if bottoming out
        print '    Prob:', n.exp((prev['score'] - step['score'])/EXP_NOISE)
        if n.random.uniform() < n.exp((prev['score'] - step['score'])/EXP_NOISE):
            # Clear off unnecessary data
            steps.append(step)
            resdat = _resdat
            if PLOT:
                _rd,_gv = 0, 0
                for bl in __resdat:
                    _rd += __resdat[bl]
                    _gv += goodval[bl]
                rd = _rd / _gv.clip(1,n.Inf)
                gv = _gv.clip(0,1)
                spec_avg = n.sum(msrdat[bl], axis=0) / n.sum(msrval[bl], axis=0).clip(1,n.Inf)
                plot_d1 = n.log10(n.abs(rd)).clip(-10,n.Inf)
                #plot_d1 = (rd*gv).real[:,590:610].flatten()
                #print n.average(plot_d1), n.std(plot_d1), plot_d1.size
                #h,bins = n.histogram(plot_d1.compress(n.abs(plot_d1) > 0), bins=100)
                #bins = (bins[1:] + bins[:-1])/2
                _d = ifft2(rd)
                #_d, info = a.deconv.clean(_d, k, tol=1e-3)
                #_d += info['res'] #/ g
                __d = a.img.recenter(_d, n.array(_d.shape)/2)
                plot_d2 = n.log10(n.abs(__d)).clip(-10,n.Inf)
                plot_d3a = n.abs(__d[gv.shape[0]/2])
                plot_d3b = n.abs(__d[gv.shape[0]/2+1])
                plot_d3c = n.abs(__d[gv.shape[0]/2-1])
                #_d *= filter
                #_d[0] = 0
                #plot_d4a = n.abs(n.average(fft2(_d), axis=0))
                #plot_d4a = n.abs(n.sum(resdat*msrval, axis=0) / n.sum(msrval, axis=0))
                #plot_d4a = n.abs(n.sum((resdat-fft2(_d))*msrval*goodval, axis=0) / n.sum(msrval*goodval, axis=0))
                #plot_d4a.shape = (plot_d4a.size/8, 8)
                #plot_d4a = n.average(plot_d4a, axis=1)
                plot_d4b = n.abs(n.sum(_rd, axis=0) / n.sum(_gv, axis=0).clip(1,n.Inf))
                plot_d4b = n.convolve(plot_d4b, n.ones(8)/8, mode='same') #/ window[0].clip(1e-4,n.Inf)
                #plot_d4b = n.abs(n.average(from_coeffs(step['coeffs']), axis=0))
                #plot_d4b.shape = (plot_d4b.size/8, 8)
                #plot_d4b = n.average(plot_d4b, axis=1)
                plot_d4c = n.abs(spec_avg) #/ window[0].clip(1e-4,n.Inf)
                #plot_d4d = plot_d4b / plot_d4c.clip(1e-4,n.Inf) * tsys
                plot_d4d = plot_d4b * 1e3
                #plot_d4c.shape = (plot_d4c.size/8, 8)
                #plot_d4c = plot_d4c[:,0]
                if LINE1 is None:
                    P.subplot(2,2,1)
                    #mx = plot_d1.max()
                    mx = -2
                    #mx = -1
                    LINE1 = P.imshow(plot_d1, vmax=mx, vmin=mx-2, aspect='auto', interpolation='nearest')
                    #LINE1 = P.imshow(plot_d1, vmax=1e-3, vmin=-1e-3, aspect='auto', interpolation='nearest')
                    P.colorbar(shrink=.5)
                    #LINE1, = P.plot(bins, h)
                    P.subplot(2,2,2)
                    #mx = plot_d2.max()
                    mx = -4
                    #mx = -3
                    LINE2 = P.imshow(plot_d2, vmax=mx, vmin=mx-2, aspect='auto', interpolation='nearest')
                    P.colorbar(shrink=.5)
                    P.subplot(2,2,3)
                    LINES[0], = P.semilogy(n.abs(plot_d3a))
                    LINES[1], = P.semilogy(n.abs(plot_d3b))
                    LINES[2], = P.semilogy(n.abs(plot_d3c))
                    P.subplot(2,2,4)
                    #LINES[10], = P.semilogy(n.abs(plot_d4a))
                    P.semilogy(tsys)
                    LINES[11], = P.semilogy(n.abs(plot_d4b))
                    LINES[12], = P.semilogy(n.abs(plot_d4c))
                    LINES[13], = P.semilogy(n.abs(plot_d4d))
                    P.ylim(1e-5,1e3)
                    P.grid()
                else:
                    LINE1.set_data(plot_d1)
                    #P.subplot(2,2,1); P.xlim(bins[0], bins[-1])
                    #LINE1.set_xdata(bins)
                    #LINE1.set_ydata(h)
                    LINE2.set_data(plot_d2)
                    LINES[0].set_ydata(plot_d3a)
                    LINES[1].set_ydata(plot_d3b)
                    LINES[2].set_ydata(plot_d3c)
                    #LINES[10].set_ydata(plot_d4a)
                    LINES[11].set_ydata(plot_d4b)
                    LINES[12].set_ydata(plot_d4c)
                    LINES[13].set_ydata(plot_d4d)
        else:
            print '    Rejecting: Score %f (%5.2f%%)' % (step['score'], n.round(100*step['score']/steps[0]['score'],2))
            # Clear off unnecessary data
            step = steps[-1].copy()
            steps.append(step)

        # Draw plots
        if PLOT: P.draw()
        print '    Score: %f (%5.2f%%)' % (steps[-1]['score'], n.round(100*steps[-1]['score']/steps[0]['score'],2))
    print 'Final Score: %f (%5.2f%%)' % (steps[-1]['score'], n.round(100*steps[-1]['score']/steps[0]['score'],2))

    coeffs = {}
    autos = {}
    for bl in steps[-1]['coeffs']:
        c[str(bl)] = steps[-1]['coeffs'][bl]
        autos[str(bl)] = __resdat[bl]
    # Pare down data and save to file
    n.savez('%s__autoinfo.npz' % (filename), times=n.array(times), chans=resdat[bl].shape[1])
    n.savez( '%s__autos.npz' % (filename), **autos)
    n.savez( '%s__coeffs.npz' % (filename), **coeffs)

