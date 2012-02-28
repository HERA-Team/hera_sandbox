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
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
o.add_option('--maxiter', dest='maxiter', type='int', default=100,
    help='Maximum number of clean iterations to allow.')
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

#FINAL_MODE = 1
FINAL_MODE = 0
EXP_NOISE = 1.0
TRIM = True

# Generate a list of files for padding
def globize(s):
    if s.isdigit(): return '*'
    return s

def gen_ddr_mask(shape, dw, drw, ratio=.25):
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

def gen_window(shape):
    w = 4 * (.5 - n.abs(.5 - n.arange(shape[0], dtype=n.float)/shape[0])).clip(0,.25)
    #P.plot(w)
    #P.show()
    #import time; time.sleep(10)
    w.shape = (w.size, 1)
    return w

def set_plot(data):
    global LINE1, LINE2
    if not PLOT: return
    if True: data = n.fft.fft(data, axis=1)
    #if False: data = a.img.recenter(data, n.array(data.shape)/2)
    if True: data = a.img.recenter(data, (0,data.shape[1]/2))
    dabs = n.log10(n.abs(data))
    dang = n.angle(data)
    if LINE1 is None:
        #P.subplot(211)
        LINE1 = P.imshow(dabs, vmax=5, vmin=1, aspect='auto')
        #P.subplot(212)
        #LINE2 = P.imshow(dang, vmax=n.pi, vmin=-n.pi, aspect='auto')
    else:
        LINE1.set_data(dabs)
        #LINE2.set_data(dang)
    P.draw()

filter, filter_cache, valid, invfilter = None, None, None, None
window, window_spec = None, None
LINE1, LINE2 = None, {}
try:
    import fftw3
    print 'Using FFTW FFT'
    fft_ddr_dat, fft_ddr_fwd, fft_ddr_rev = None, None, None
    def ddr_filter(d, dw, drw):
        global LINE1, filter_cache, filter, invfilter, valid
        global window, window_spec
        global fft_ddr_dat, fft_ddr_fwd, fft_ddr_rev
        if fft_ddr_dat is None or d.shape != fft_ddr_dat.shape:
            print 'Reseting FFTW Plan'
            fft_ddr_dat = n.zeros((d.shape[0],d.shape[1]), dtype=n.complex)
            fft_ddr_fwd = fftw3.Plan(fft_ddr_dat, None, direction='forward', flags=['measure'])
            fft_ddr_rev = fftw3.Plan(fft_ddr_dat, None, direction='backward', flags=['measure'])
            window = gen_window(d.shape)
            window_spec = n.fft.fft(window, axis=0) / window.size
        #fft_ddr_dat[:] = d
        fft_ddr_dat[:] = d * window
        fft_ddr_rev() # DDR Transform
        #fft_ddr_dat /= window_spec
        #set_plot(fft_ddr_dat)
        # Apply DDR Filter
        if (dw,drw) != filter_cache:
            filter_cache = (dw, drw)
            filter = gen_ddr_mask(fft_ddr_dat.shape, dw, drw)
            valid = n.where(filter)
            invfilter = 1-filter
        if False: # Add noise for true MCMC
            #sig = n.sqrt(n.sum(n.abs(fft_ddr_dat*invfilter)[1:,1:]**2)/n.sum(invfilter[1:,1:])),
            sig = n.sqrt(fft_ddr_dat.size) * EXP_NOISE
            #print '    Sigma:', sig
            noise_amp = n.random.normal(scale=sig, size=valid[0].shape[0])
            noise_phs = n.random.uniform(0,2*n.pi, size=valid[0].shape[0])
            noise = noise_amp * n.exp(1j*noise_phs)
            fft_ddr_dat[valid] += noise
        #fft_ddr_dat[:] = fft_ddr_dat * filter
        print n.abs(fft_ddr_dat[valid])
        fft_ddr_dat *= filter
        #print n.sqrt(n.sum(n.abs(fft_ddr_dat)**2)/n.sum(filter))
        fft_ddr_fwd() # undo DDR Transform
        return fft_ddr_dat / fft_ddr_dat.size
except(ImportError):
    print 'Using numpy FFT'
    def ddr_filter(d, dw, drw):
        global LINE1, filter_cache, filter
        global window, window_spec
        # Apply DDR Filter
        if (dw,drw) != filter_cache:
            filter_cache = (dw, drw)
            filter = gen_ddr_mask(d.shape, dw, drw)
            window = gen_window(d.shape)
            window_spec = n.fft.fft(window, axis=0) / window.size
        d = n.fft.ifft2(d * window) # DDR Transform
        d *= filter
        d = n.fft.fft2(d) # undo DDR Transform
        return d

def compute_residual(msrdat, msrval, simdat, rsvdat, __srcest_ant, __srcest_bl, xtalk, do_xtalk=True):
    # Compute residual for all sources
    __resdat = {}
    for bl in msrdat:
        i,j = a.miriad.bl2ij(bl)
        __resdat[bl] = msrdat[bl].copy()
        for k in simdat:
            if rsvdat.has_key(k): mdl = simdat[k][bl] * rsvdat[k][bl]
            else: mdl = simdat[k][bl]
            __resdat[bl] -= __srcest_ant.get(k,{}).get(i,0) * n.conj(__srcest_ant.get(k,{}).get(j,0)) * mdl
            if __srcest_bl.has_key(k): __resdat[bl] -= __srcest_bl[k][bl] * simdat[k][bl]
        if not do_xtalk: __resdat[bl] -= xtalk[bl]
        __resdat[bl] *= msrval[bl] # Mask out invalid data
        #if bl == a.miriad.ij2bl(0,3):
        #if bl == a.miriad.ij2bl(0,14):
            #set_plot(__resdat[bl] * window)
    if do_xtalk: return compute_xtalk(msrval, __resdat)
    else: return __resdat

def compute_xtalk(msrval, __resdat):
    __xtalk = {}
    for bl in msrval:
        xsum = n.sum(__resdat[bl], axis=0)
        xwgt = n.sum(msrval[bl], axis=0)
        __xtalk[bl] = xsum/xwgt.clip(1,n.Inf)
        __resdat[bl] -= __xtalk[bl]
        __resdat[bl] *= msrval[bl] # Mask out invalid data
    return __resdat, __xtalk

def compute_score(resdat, msrval):
    return n.sqrt(sum([n.sum(n.abs((resdat[bl]*window)**2)) for bl in msrval]) / (sum([n.sum(msrval[bl]*window) for bl in msrval])))

file_glob = '.'.join(map(globize, args[0].split('.')))
filelist = glob.glob(file_glob)
filelist.sort()

for arg in args:
    # Gather data
    print 'Processing', arg
    findex = filelist.index(arg)
    files = filelist[findex:findex+1]
    #files = filelist[findex:findex+2]
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
    iscore = score = n.sqrt(sum([n.sum(n.abs(msrdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))
    print 'Initial Score: %f' % (score)
    hist = {}
    #dw,drw = 0,0
    dw,drw = opts.dw, opts.drw
    mode,tier = -1,0
    srcest_ant, srcest_bl = {}, {}
    xtalk = {}
    # Create residuals
    resdat = {}
    for bl in msrdat: resdat[bl] = msrdat[bl].copy()

    for iter in xrange(opts.maxiter):
        print 'Iteration %d:    DW=%d, DRW=%d, MODE=%d, TIER=%d' % (iter, dw, drw, mode, tier)
        # Derive source estimators
        _srcest_ant, _srcest_bl = {}, {}
        if mode == -1: # BEAMFORM
            for k in simdat:
                if not k in srctier[tier]:
                    try: _srcest_ant[k] = srcest_ant[k]
                    except(KeyError): pass
                    continue
                d, w = 0, 0
                for bl in simdat[k]:
                    if k in blsrcs: wgt = blwgt[k][bl] * rsvdat[k][bl]
                    else: wgt = blwgt[k][bl]
                    d += (resdat[bl] * n.conj(simdat[k][bl])).real * wgt
                    #d += (resdat[bl] * n.conj(simdat[k][bl])) * wgt
                    if k in blsrcs: wgt *= rsvdat[k][bl]
                    w += wgt
                d /= n.where(w == 0, 1, w)
                d = ddr_filter(n.where(w == 0, 0, n.sqrt(n.abs(d))), dw, drw)
                #d = ddr_filter(d, dw, drw)
                #d = n.where(w == 0, 0, n.sqrt(d))
                if not _srcest_ant.has_key(k): _srcest_ant[k] = {}
                for ant in ants:
                    if not srcest_ant.get(k,{}).has_key(ant): _srcest_ant[k][ant] = d
                    else: _srcest_ant[k][ant] = srcest_ant[k][ant]
            mode = 0
        elif mode == 0: # VARY PER ANT
            for k in simdat:
                if not srcest_ant.has_key(k): continue
                if not k in srctier[tier]:
                    try: _srcest_ant[k] = srcest_ant[k]
                    except(KeyError): pass
                    continue
                d,w = {}, {}
                for bl in simdat[k]:
                    i,j = a.miriad.bl2ij(bl)
                    if k in blsrcs: wgt = blwgt[k][bl] * rsvdat[k][bl]
                    else: wgt = blwgt[k][bl]
                    _d = resdat[bl] * n.conj(simdat[k][bl]) * wgt
                    if k in blsrcs: wgt *= rsvdat[k][bl]
                    _w = wgt
                    d[i], w[i] = d.get(i,0) + .5*_d, w.get(i,0) + _w * n.conj(srcest_ant.get(k,{}).get(j,0))
                    d[j], w[j] = d.get(j,0) + n.conj(.5*_d), w.get(j,0) + _w * n.conj(srcest_ant.get(k,{}).get(i,0))
                for i in d:
                    # This divide can cause an explosion if w[i] is close to
                    # 0, which can happen early on when summing weights that
                    # depend on estimated antenna gain.  Adding noise helps
                    d[i] /= n.where(w[i] == 0, 1, w[i])
                    d[i] = ddr_filter(d[i], dw, drw)
                    d[i] = n.where(w[i] == 0, 0, d[i])
                    if not srcest_ant.has_key(k): srcest_ant[k] = {}
                    if not _srcest_ant.has_key(k): _srcest_ant[k] = {}
                    #if i == 0 and k == '0:42:53.96_52:07:34.8':
                    #    P.subplot(234)
                    #    P.imshow(n.log10(n.abs(d[i])), aspect='auto')
                    #    #P.imshow(a.img.recenter(n.log10(n.abs(n.fft.fft2(d[i]))), n.array(d[i].shape)/2), aspect='auto')
                    #    #P.imshow(n.angle(d[i]), aspect='auto')
                    _srcest_ant[k][i] = srcest_ant[k].get(i,0) + d[i]
                    #if i == 0 and k == '0:42:53.96_52:07:34.8':
                    #    P.subplot(235)
                    #    P.imshow(n.log10(n.abs(srcest_ant[k][i])), aspect='auto')
                    #    #P.imshow(a.img.recenter(n.log10(n.abs(n.fft.fft2(srcest_ant[k][i]))), n.array(d[i].shape)/2), aspect='auto')
                    #    #P.imshow(n.angle(srcest_ant[k][i]), aspect='auto')
                    #    P.subplot(236)
                    #    P.imshow(n.log10(n.abs(_srcest_ant[k][i])), aspect='auto')
                    #    #P.imshow(a.img.recenter(n.log10(n.abs(n.fft.fft2(_srcest_ant[k][i]))), n.array(d[i].shape)/2), aspect='auto')
                    #    #P.imshow(n.angle(_srcest_ant[k][i]), aspect='auto')
                    #    #P.colorbar()
                    #    #import time; time.sleep(10)
        else: # Polish up per baseline for select sources
            for k in simdat:
                if not k in blsrcs: continue
                if not k in srctier[tier]:
                    try: _srcest_bl[k] = srcest_bl[k]
                    except(KeyError): pass
                    continue
                for bl in msrdat:
                    d = resdat[bl] * n.conj(simdat[k][bl])
                    d = ddr_filter(d, dw, drw)
                    d = n.where(blwgt[k][bl] == 0, 0, d)
                    if not _srcest_bl.has_key(k): _srcest_bl[k] = {}
                    _srcest_bl[k][bl] = srcest_bl.get(k,{}).get(bl,0) + d

        # Compute residual for all sources
        if mode == 0:
            _resdat,_xtalk = compute_residual(msrdat, msrval, simdat, rsvdat, _srcest_ant, srcest_bl, xtalk, do_xtalk=True)
        else:
            _resdat = compute_residual(msrdat, msrval, simdat, rsvdat, srcest_ant, _srcest_bl, xtalk, do_xtalk=False)
            _xtalk = xtalk
        _score = compute_score(_resdat,msrval)

        # Always accept any downhill jump, advance state if bottoming out
        tol = n.abs(1 - _score/score)
        #tol = 1 - _score/score
        if n.random.uniform() < n.exp((score - _score)/EXP_NOISE):
            resdat, score, xtalk = _resdat, _score, _xtalk
            if mode <= 0: srcest_ant = _srcest_ant
            else: srcest_bl = _srcest_bl
            h, bins = n.histogram([score], bins=1000, range=(0.,20.))
            hist['cnt'] = hist.get('cnt',0) + 1.
            hist['score'] = hist.get('score',0) + h.astype(n.float)
            for k in srcest_ant:
                d = n.average(n.abs(srcest_ant[k][0])**2, axis=1)
                #hist[k] = hist.get(k,0) + d
                hist[k] = d
                #if k == '0:42:53.96_52:07:34.8':
                #    P.subplot(133)
                #    P.imshow(n.abs(srcest_ant[k][0])**2, aspect='auto')
            if LINE1 is None:
                P.subplot(121)
                bins = .5*(bins[1:] + bins[:-1])
                LINE1, = P.plot(bins, hist['score']/hist['cnt'], 'k-')
                P.xlim(0,20)
                P.ylim(0,1)
                P.subplot(122)
                for k in srcest_ant:
                    #L, = P.semilogy(hist[k]/hist['cnt'], '-')
                    #if n.average(hist[k]) > 1: continue
                    if TRIM:
                        i0 = len(times) / 6
                        LINE2[k], = P.semilogy(hist[k][i0:-i0], '-', label=k)
                    else:
                        LINE2[k], = P.semilogy(hist[k], '-', label=k)
                    #print k, n.sum(n.abs((d[1:-1] - .5*(d[2:]+d[:-2]))/d[1:-1])**2)
                P.ylim(.1,1e4)
                #P.legend()
            else:
                LINE1.set_ydata(hist['score']/hist['cnt'])
                for k in LINE2:
                    if TRIM:
                        i0 = len(times) / 6
                        LINE2[k].set_ydata(hist[k][i0:-i0])
                    else:
                        #LINE2[k].set_ydata(hist[k]/hist['cnt'])
                        LINE2[k].set_ydata(hist[k])
        else:
            print '    Rejecting proposal: Score %f (%5.2f%%, tol=%e)' % (_score, n.round(100*_score/iscore,2), tol)

        P.draw()
        # If we are at tolerance, add degrees of freedom and repeat
        if tol < opts.clean:
            if mode < FINAL_MODE:
                mode += 1
            #else:
            #    if dw == opts.dw and drw == opts.drw:
            #        if tier+1 < len(srctier):
            #            tier += 1
            #            #mode = 0
            #            mode = -1
            #            dw,drw = 0,0
            #        else:
            #            print '    Tolerance'
            #            break
            #    else:
            #        mode = 0
            #        dw += 1; drw += 1
        print '    New Score: %f (%5.2f%%, tol=%e)' % (score, n.round(100*score/iscore,2), tol)
    print '    Final Score: %f (%5.2f%%, tol=%f)' % (score, n.round(100*score/iscore,2), tol)

    # Pare down data and save to file
    if TRIM:
        i0 = len(times) / 6
        times = times[i0:-i0]
    n.savez('%s__times.npz' % (arg), times=n.array(times))
    n.savez('%s__afreqs.npz' % (arg), freqs=afreqs)
    __xtalk = {}
    for bl in xtalk: __xtalk[str(bl)] = xtalk[bl]
    n.savez( '%s__xtalk.npz' % (arg), **__xtalk)
    if FINAL_MODE >= 0:
        for k in srcest_ant:
            d = {}
            for i in srcest_ant[k]:
                if TRIM: d[str(i)] = srcest_ant[k][i][i0:-i0]
                else: d[str(i)] = srcest_ant[k][i]
            n.savez( '%s__srcest_ant__%s.npz' % (arg,k), **d)
    if FINAL_MODE >= 1:
        for k in srcest_bl:
            d = {}
            for bl in srcest_bl[k]:
                if TRIM: d[str(bl)] = srcest_bl[k][bl][i0:-i0]
                else: d[str(bl)] = srcest_bl[k][bl]
            n.savez( '%s__srcest_bl__%s.npz' % (arg,k), **d)

