#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, math, glob
PLOT = False
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

CLEAN_GAIN = .9
#CLEAN_GAIN = .5
XTALK_GAIN = 1
FINAL_MODE = 1
#FINAL_MODE = 0

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

def set_plot(data):
    global LINE1
    if not PLOT: return
    if LINE1 is None:
        LINE1 = P.imshow(a.img.recenter(n.log10(n.abs(data)), n.array(data.shape)/2), vmax=7, vmin=0)
    else:
        LINE1.set_data(a.img.recenter(n.log10(n.abs(data)), n.array(data.shape)/2))
    P.draw()

filter, filter_cache = None, None
LINE1 = None
try:
    import fftw3
    print 'Using FFTW FFT'
    fft_ddr_dat, fft_ddr_fwd, fft_ddr_rev = None, None, None
    def ddr_filter(d, dw, drw):
        global LINE1, filter_cache, filter
        global fft_ddr_dat, fft_ddr_fwd, fft_ddr_rev
        #padlen = math.ceil(d.shape[1] * .2)
        #d = n.concatenate([n.fliplr(d[:,:padlen]),d,n.fliplr(d[:,-padlen:])], axis=1)
        if fft_ddr_dat is None or d.shape != fft_ddr_dat.shape:
            print 'Reseting FFTW Plan'
            fft_ddr_dat = n.zeros((d.shape[0],d.shape[1]), dtype=n.complex)
            fft_ddr_fwd = fftw3.Plan(fft_ddr_dat, None, direction='forward', flags=['measure'])
            fft_ddr_rev = fftw3.Plan(fft_ddr_dat, None, direction='backward', flags=['measure'])
        fft_ddr_dat[:] = d
        fft_ddr_rev() # DDR Transform
        set_plot(fft_ddr_dat)
        # Apply DDR Filter
        if (dw,drw) != filter_cache:
            filter_cache = (dw, drw)
            filter = gen_ddr_mask(fft_ddr_dat.shape, dw, drw)
        fft_ddr_dat[:] = fft_ddr_dat * filter
        fft_ddr_fwd() # undo DDR Transform
        #return fft_ddr_dat[:,padlen:-padlen] / (d.shape[0] * d.shape[1])
        return fft_ddr_dat / fft_ddr_dat.size #(d.shape[0] * d.shape[1])
except(ImportError):
    print 'Using numpy FFT'
    def ddr_filter(d, dw, drw):
        global LINE1, filter_cache, filter
        padlen = math.ceil(d.shape[1] * .2)
        d = n.fft.ifft2(d) # DDR Transform
        # Apply DDR Filter
        if (dw,drw) != filter_cache:
            filter_cache = (dw, drw)
            filter = gen_ddr_mask(fft_ddr_dat.shape, dw, drw)
        d *= filter
        d = n.fft.fft2(d) # undo DDR Transform
        return d#[:,padlen:-padlen]


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
                    #simd = n.conj(aa.gen_phs(cat[k],i,j, srcshape=cat[k].srcshape, resolve_src=True))
                    simd = n.conj(aa.gen_phs(cat[k],i,j))
                    u,v,w = aa.gen_uvw(i,j,cat[k]).squeeze()
                    if k in blsrcs: rsvd = aa.resolve_src(u, v, srcshape=cat[k].srcshape)
                except(a.phs.PointingError):
                    simd = n.zeros_like(d)
                    u = v = w = n.zeros_like(d)
                    if k in blsrcs: rsvd = simd
                simdat[k][bl] = simdat[k].get(bl,[]) + [simd]
                if k in blsrcs: rsvdat[k][bl] = rsvdat[k].get(bl,[]) + [rsvd]
                blwgt[k][bl] = f * n.sqrt(u**2 + v**2)
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
    dw,drw = 0,0
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
                    if k in blsrcs: wgt *= rsvdat[k][bl]
                    w += wgt
                d /= n.where(w > 0, w, 1)
                d = ddr_filter(d, dw, drw)
                d = n.sqrt(d)
                if not _srcest_ant.has_key(k): _srcest_ant[k] = {}
                for ant in ants:
                    if not srcest_ant.get(k,{}).has_key(ant): _srcest_ant[k][ant] = CLEAN_GAIN * d
                    else: _srcest_ant[k][ant] = srcest_ant[k][ant]
                    #_srcest_ant[k][ant] = srcest_ant.get(k,{}).get(ant,0) + CLEAN_GAIN * d / srcest_ant.get(k,{}).get(ant,n.sqrt(d))
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
                    #d[i], w[i] = d.get(i,0) + _d, w.get(i,0) + _w*(est_gi+n.conj(est_gj))
                    #d[j], w[j] = d.get(j,0) + n.conj(_d), w.get(j,0) + _w*(n.conj(est_gi)+est_gj)
                    d[i], w[i] = d.get(i,0) + .5*_d, w.get(i,0) + _w * n.conj(srcest_ant.get(k,{}).get(j,0))
                    d[j], w[j] = d.get(j,0) + n.conj(.5*_d), w.get(j,0) + _w * n.conj(srcest_ant.get(k,{}).get(i,0))
                for i in d:
                    d[i] /= n.where(w[i] == 0, 1, w[i])
                    d[i] = ddr_filter(d[i], dw, drw)
                    if not srcest_ant.has_key(k): srcest_ant[k] = {}
                    if not _srcest_ant.has_key(k): _srcest_ant[k] = {}
                    _srcest_ant[k][i] = srcest_ant[k].get(i,0) + CLEAN_GAIN * d[i]
        else: # Polish up per baseline for select sources
            for k in simdat:
                if not k in blsrcs: continue
                if not k in srctier[tier]:
                    try: _srcest_bl[k] = srcest_bl[k]
                    except(KeyError): pass
                    continue
                for bl in msrdat:
                    #sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                    #d = resdat[bl] * n.conj(simdat[k][bl]) / sd
                    d = resdat[bl] * n.conj(simdat[k][bl])
                    #if bl == a.miriad.ij2bl(0,8): PLOT = True
                    #else: PLOT = False
                    d = ddr_filter(d, dw, drw)
                    #if not srcest_bl.has_key(k): srcest_bl[k] = {}
                    if not _srcest_bl.has_key(k): _srcest_bl[k] = {}
                    _srcest_bl[k][bl] = srcest_bl.get(k,{}).get(bl,0) + CLEAN_GAIN * d
                    #print k, a.miriad.bl2ij(bl), n.sqrt(n.average(n.abs(_srcest_bl[k][bl])**2))
        # Estimate XTALK too
        _xtalk = {}
        for bl in msrdat:
            xsum = n.sum(resdat[bl], axis=0)
            xwgt = n.sum(msrval[bl], axis=0)
            _xtalk[bl] = xtalk.get(bl,0) + XTALK_GAIN * xsum/xwgt.clip(1,n.Inf)

        # Compute residual for all sources
        _resdat = {}
        for bl in msrdat:
            i,j = a.miriad.bl2ij(bl)
            _resdat[bl] = msrdat[bl].copy()
            for k in simdat:
                if k in blsrcs: mdl = simdat[k][bl] * rsvdat[k][bl]
                else: mdl = simdat[k][bl]
                if mode == 0: _resdat[bl] -= _srcest_ant.get(k,{}).get(i,0) * n.conj(_srcest_ant.get(k,{}).get(j,0)) * mdl
                else: _resdat[bl] -= srcest_ant.get(k,{}).get(i,0) * n.conj(srcest_ant.get(k,{}).get(j,0)) * mdl
                if mode == 1 and _srcest_bl.has_key(k):
                    #phsdat = simdat[k][bl] / n.abs(simdat[k][bl])
                    #sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                    #_resdat[bl] -= _srcest_bl[k][bl] * simdat[k][bl] / sd
                    _resdat[bl] -= _srcest_bl[k][bl] * simdat[k][bl]
                elif srcest_bl.has_key(k):
                    #phsdat = simdat[k][bl] / n.abs(simdat[k][bl])
                    #sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                    #_resdat[bl] -= srcest_bl[k][bl] * simdat[k][bl] /sd
                    _resdat[bl] -= srcest_bl[k][bl] * simdat[k][bl]
            _resdat[bl] -= _xtalk[bl]
            _resdat[bl] *= msrval[bl] # Mask out invalid data
        _score = n.sqrt(sum([n.sum(n.abs(_resdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))

        # Recompute source by source if this score is not an improvement
        if _score > score:
            # Figure out which changes to model improve the residuals
            for k in simdat:
                if k not in srctier[tier]: continue
                if mode == 1 and not k in blsrcs: continue
                __resdat = {}
                for bl in msrdat:
                    i,j = a.miriad.bl2ij(bl)
                    __resdat[bl] = resdat[bl].copy()
                    if k in blsrcs: mdl = simdat[k][bl] * rsvdat[k][bl]
                    else: mdl = simdat[k][bl]
                    __resdat[bl] -= (_srcest_ant[k][i]*n.conj(_srcest_ant[k][j])-srcest_ant[k][i]*n.conj(srcest_ant[k][j])) * mdl
                    if _srcest_bl.has_key(k):
                        #phsdat = simdat[k][bl] / n.abs(simdat[k][bl])
                        #sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                        #__resdat[bl] -= (_srcest_bl[k][bl] - srcest_bl.get(k,{}).get(bl,0)) * simdat[k][bl] / sd
                        __resdat[bl] -= (_srcest_bl[k][bl] - srcest_bl.get(k,{}).get(bl,0)) * simdat[k][bl]
                    __resdat[bl] -= (_xtalk[bl] - xtalk.get(bl,0))
                    __resdat[bl] *= msrval[bl] # Mask out invalid data
                __score = n.sqrt(sum([n.sum(n.abs(__resdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))
                if __score >= score: 
                    print '      * %16s %f' % (k, __score-score)
                    if mode <= 0:
                        #if srcest_ant.has_key(k): _srcest_ant[k] = srcest_ant[k].copy()
                        if srcest_ant.has_key(k): _srcest_ant[k] = srcest_ant[k]
                        elif _srcest_ant.has_key(k): del(_srcest_ant[k])
                    else:
                        for bl in _srcest_bl[k]: _srcest_bl[k][bl] = srcest_bl.get(k,{}).get(bl,0)
                else: print '    %20s %f' % (k, __score-score)

            # Compute residual for sources that succeeded individually
            _resdat = {}
            for bl in msrdat:
                i,j = a.miriad.bl2ij(bl)
                _resdat[bl] = msrdat[bl].copy()
                for k in simdat:
                    if k in blsrcs: mdl = simdat[k][bl] * rsvdat[k][bl]
                    else: mdl = simdat[k][bl]
                    _resdat[bl] -= _srcest_ant.get(k,{}).get(i,0) * n.conj(_srcest_ant.get(k,{}).get(j,0)) * mdl
                    if mode == 1 and _srcest_bl.has_key(k):
                        #phsdat = simdat[k][bl] / n.abs(simdat[k][bl])
                        #sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                        #_resdat[bl] -= _srcest_bl[k][bl] * simdat[k][bl] / sd
                        _resdat[bl] -= _srcest_bl[k][bl] * simdat[k][bl]
                    elif srcest_bl.has_key(k):
                        #phsdat = simdat[k][bl] / n.abs(simdat[k][bl])
                        #sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                        #_resdat[bl] -= srcest_bl[k][bl] * simdat[k][bl] /sd
                        _resdat[bl] -= srcest_bl[k][bl] * simdat[k][bl]
                _resdat[bl] -= _xtalk[bl]
                _resdat[bl] *= msrval[bl] # Mask out invalid data
            _score = n.sqrt(sum([n.sum(n.abs(_resdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))

        # Always accept any downhill jump, advance state if bottoming out
        #tol = 1 - _score/score
        tol = n.abs(1 - _score/score)
        if tol < 0:
            print '    Divergence'
            print '    Failed Score: %f (%5.2f%%)' % (_score, n.round(100*_score/iscore,2))
        else: 
            resdat, score, xtalk = _resdat, _score, _xtalk
            if mode <= 0: srcest_ant = _srcest_ant
            else: srcest_bl = _srcest_bl

        # If we are at tolerance, add degrees of freedom and repeat
        if tol < opts.clean:
            if mode < FINAL_MODE:
                mode += 1
            else:
                if dw == opts.dw and drw == opts.drw:
                    if tier+1 < len(srctier):
                        tier += 1
                        #mode = 0
                        mode = -1
                        dw,drw = 0,0
                    else:
                        print '    Tolerance'
                        break
                else:
                    mode = 0
                    dw += 1; drw += 1
        print '    New Score: %f (%5.2f%%, tol=%e)' % (score, n.round(100*score/iscore,2), tol)
    print '    Final Score: %f (%5.2f%%, tol=%f)' % (score, n.round(100*score/iscore,2), tol)

    # Pare down data and save to file
    #i0 = len(times) / 6
    #times = times[i0:-i0]
    n.savez('%s__times.npz' % (arg), times=n.array(times))
    n.savez('%s__afreqs.npz' % (arg), freqs=afreqs)
    __xtalk = {}
    for bl in xtalk: __xtalk[str(bl)] = xtalk[bl]
    n.savez( '%s__xtalk.npz' % (arg), **__xtalk)
    if FINAL_MODE >= 0:
        for k in srcest_ant:
            d = {}
            for i in srcest_ant[k]:
                #d[str(i)] = srcest_ant[k][i][i0:-i0]
                d[str(i)] = srcest_ant[k][i]
            n.savez( '%s__srcest_ant__%s.npz' % (arg,k), **d)
    if FINAL_MODE >= 1:
        for k in srcest_bl:
            d = {}
            for bl in srcest_bl[k]:
                #d[str(bl)] = srcest_bl[k][bl][i0:-i0]
                d[str(bl)] = srcest_bl[k][bl]
            n.savez( '%s__srcest_bl__%s.npz' % (arg,k), **d)

