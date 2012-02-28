#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, math, glob

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

CLEAN_GAIN = .2
XTALK_GAIN = 1
FINAL_MODE = 2

# Generate a list of files for padding
def globize(s):
    if s.isdigit(): return '*'
    return s

def ddr_filter(d, drw, dw):
    padlen = math.ceil(d.shape[1] * .2)
    d = n.concatenate([n.fliplr(d[:,:padlen]),d,n.fliplr(d[:,-padlen:])], axis=1)
    d = n.fft.ifft(d, axis=1) # Delay Transform
    d = n.fft.ifft(d, axis=0) # Delay-Rate Transform
    # Apply DDR Filter
    x1,x2 = drw, -drw
    if x2 == 0: x2 = d.shape[0]
    y1,y2 = dw, -dw
    if y2 == 0: y2 = d.shape[1]
    d[x1+1:x2] = 0
    d[:,y1+1:y2] = 0
    d = n.fft.fft(d, axis=0) # undo Delay-Rate Transform
    d = n.fft.fft(d, axis=1) # undo Delay Transform
    return d[:,padlen:-padlen]

file_glob = '.'.join(map(globize, args[0].split('.')))
filelist = glob.glob(file_glob)
filelist.sort()

for arg in args:
    # Gather data
    print 'Processing', arg
    findex = filelist.index(arg)
    #files = filelist[findex:findex+1]
    #files = filelist[findex:findex+2]
    files = filelist[findex-1:findex+2]
    msrdat, msrval, simdat = {}, {}, {}
    blwgt = {}
    for k in cat:
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
            for k in simdat.keys():
                try:
                    simd = n.conj(aa.gen_phs(cat[k],i,j, 
                        srcshape=cat[k].srcshape, resolve_src=True))
                    u,v,w = aa.gen_uvw(i,j,cat[k]).squeeze()
                except(a.phs.PointingError):
                    simd = n.zeros_like(d)
                    u = v = w = n.zeros_like(d)
                simdat[k][bl] = simdat[k].get(bl,[]) + [simd]
                blwgt[k][bl] = f * n.sqrt(u**2 + v**2)
    # Simulate visibilities for each source for the data we collected
    for bl in msrdat:
        msrdat[bl] = n.array(msrdat[bl])
        msrval[bl] = n.array(msrval[bl])
        for k in simdat.keys():
            simdat[k][bl] = n.array(simdat[k][bl])
            blwgt[k][bl] = n.array(blwgt[k][bl])
            if n.all(simdat[k][bl] == 0):
                del(simdat[k])
                del(blwgt[k])
    iscore = score = n.sqrt(sum([n.sum(n.abs(msrdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))
    print 'Initial Score: %f' % (score)
    dw,drw = 0,0
    mode,tier = 0,0
    srcest_bm, srcest_ant, srcest_bl = {}, {}, {}
    xtalk = {}
    # Create residuals
    resdat = {}
    for bl in msrdat: resdat[bl] = msrdat[bl].copy()

    for iter in xrange(opts.maxiter):
        print 'Iteration %d:    DW=%d, DRW=%d, MODE=%d, TIER=%d' % (iter, dw, drw, mode, tier)
        # Derive source estimators
        _srcest_bm, _srcest_ant, _srcest_bl = {}, {}, {}
        if mode == 0: # BEAMFORM
            for k in simdat:
                if not k in srctier[tier]:
                    try: _srcest_bm[k] = srcest_bm[k]
                    except(KeyError): pass
                    continue
                d, w = 0., 0.
                for bl in simdat[k]:
                    d += resdat[bl] * n.conj(simdat[k][bl]) * blwgt[k][bl]
                    w += n.abs(simdat[k][bl])**2 * blwgt[k][bl]
                d /= n.where(w > 0, w, 1)
                d = ddr_filter(d, dw, drw)
                _srcest_bm[k] = srcest_bm.get(k, 0.) + CLEAN_GAIN * d
        elif mode == 1: # VARY PER ANT
            for k in simdat:
                if not srcest_bm.has_key(k): continue
                if not k in srctier[tier]:
                    try: _srcest_ant[k] = srcest_ant[k]
                    except(KeyError): pass
                    continue
                d,w = {}, {}
                for bl in simdat[k]:
                    i,j = a.miriad.bl2ij(bl)
                    _d = resdat[bl] * n.conj(simdat[k][bl]) * blwgt[k][bl]
                    _w = n.abs(simdat[k][bl])**2 * blwgt[k][bl]
                    est_gi = n.sqrt(srcest_bm[k]) + srcest_ant.get(k,{}).get(i,0.)
                    est_gj = n.sqrt(srcest_bm[k]) + srcest_ant.get(k,{}).get(j,0.)
                    d[i], w[i] = d.get(i,0.) + _d, w.get(i,0.) + _w*(est_gi+n.conj(est_gj))
                    d[j], w[j] = d.get(j,0.) + n.conj(_d), w.get(j,0.) + _w*(n.conj(est_gi)+est_gj)
                for i in d:
                    d[i] /= n.where(w[i] > 0, w[i], 1)
                    d[i] = ddr_filter(d[i], dw, drw)
                    if not srcest_ant.has_key(k): srcest_ant[k] = {}
                    if not _srcest_ant.has_key(k): _srcest_ant[k] = {}
                    _srcest_ant[k][i] = srcest_ant[k].get(i,0.) + CLEAN_GAIN * d[i]
        else: # Polish up per baseline for select sources
            for k in simdat:
                if not k in blsrcs: continue
                if not k in srctier[tier]:
                    try: _srcest_bl[k] = srcest_bl[k]
                    except(KeyError): pass
                    continue
                for bl in msrdat:
                    sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                    d = resdat[bl] * n.conj(simdat[k][bl]) / sd
                    d = ddr_filter(d, dw, drw)
                    if not srcest_bl.has_key(k): srcest_bl[k] = {}
                    if not _srcest_bl.has_key(k): _srcest_bl[k] = {}
                    _srcest_bl[k][bl] = srcest_bl[k].get(bl,0.) + CLEAN_GAIN * d
        # Estimate XTALK too
        _xtalk = {}
        for bl in msrdat:
            xsum = n.sum(resdat[bl], axis=0)
            xwgt = n.sum(msrval[bl], axis=0)
            _xtalk[bl] = xtalk.get(bl,0.) + XTALK_GAIN * xsum/xwgt.clip(1,n.Inf)

        # Figure out which changes to model improve the residuals
        for k in simdat:
            if k not in srctier[tier]: continue
            if mode == 2 and not k in blsrcs: continue
            __resdat = {}
            for bl in msrdat:
                i,j = a.miriad.bl2ij(bl)
                __resdat[bl] = resdat[bl].copy()
                gi,gj = 0,0
                _gi,_gj = 0,0
                if srcest_bm.has_key(k):
                    gi = n.sqrt(srcest_bm[k])
                    gj = gi.copy()
                    _gi,_gj = gi.copy(),gj.copy()
                if _srcest_bm.has_key(k):
                    _gi = n.sqrt(_srcest_bm[k])
                    _gj = _gi.copy()
                if srcest_ant.get(k,{}).has_key(i): gi += srcest_ant[k][i]
                if srcest_ant.get(k,{}).has_key(j): gj += srcest_ant[k][j]
                if _srcest_ant.get(k,{}).has_key(i): _gi += _srcest_ant[k][i]
                elif srcest_ant.get(k,{}).has_key(i): _gi += srcest_ant[k][i]
                if _srcest_ant.get(k,{}).has_key(j): _gj += _srcest_ant[k][j]
                elif srcest_ant.get(k,{}).has_key(j): _gj += srcest_ant[k][j]
                # Subtract only the difference between the proposed model
                # and the model that's been subtracted already
                __resdat[bl] -= (_gi*n.conj(_gj)-gi*n.conj(gj)) * simdat[k][bl]
                sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                if _srcest_bl.has_key(k):
                    __resdat[bl] -= (_srcest_bl[k][bl] - srcest_bl.get(k,{}).get(bl,0.)) * simdat[k][bl] / sd
                __resdat[bl] -= (_xtalk[bl] - xtalk.get(bl,0.))
                __resdat[bl] *= msrval[bl] # Mask out invalid data
            __score = n.sqrt(sum([n.sum(n.abs(__resdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))
            if __score >= score: 
                print '      * %16s %f' % (k, __score-score)
                if mode == 0:
                    if srcest_bm.has_key(k): _srcest_bm[k] = srcest_bm[k].copy()
                    elif _srcest_bm.has_key(k): del(_srcest_bm[k])
                elif mode == 1:
                    if srcest_ant.has_key(k): _srcest_ant[k] = srcest_ant[k].copy()
                    elif _srcest_ant.has_key(k): del(_srcest_ant[k])
                else:
                    for bl in _srcest_bl[k]: _srcest_bl[k][bl] = srcest_bl.get(k,{}).get(bl,0.)
            else:
                print '    %20s %f' % (k, __score-score)

        # Compute residual for sources that succeeded individually
        _resdat = {}
        for bl in msrdat:
            i,j = a.miriad.bl2ij(bl)
            _resdat[bl] = msrdat[bl].copy()
            for k in simdat:
                gi,gj = 0,0
                if mode == 0:
                    if _srcest_bm.has_key(k):
                        gi = n.sqrt(_srcest_bm[k])
                        gj = gi.copy()
                elif srcest_bm.has_key(k):
                    gi = n.sqrt(srcest_bm[k])
                    gj = gi.copy()
                if mode == 1:
                    if _srcest_ant.get(k,{}).has_key(i): gi += _srcest_ant[k][i]
                    if _srcest_ant.get(k,{}).has_key(j): gj += _srcest_ant[k][j]
                elif srcest_ant.has_key(k):
                    if srcest_ant[k].has_key(i): gi += srcest_ant[k][i]
                    if srcest_ant[k].has_key(j): gj += srcest_ant[k][j]
                _resdat[bl] -= gi * n.conj(gj) * simdat[k][bl]
                sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                if mode == 2 and _srcest_bl.has_key(k): _resdat[bl] -= _srcest_bl[k][bl] * simdat[k][bl] / sd
                elif srcest_bl.has_key(k): _resdat[bl] -= srcest_bl[k][bl] * simdat[k][bl] /sd
            _resdat[bl] -= _xtalk[bl]
            _resdat[bl] *= msrval[bl] # Mask out invalid data
        _score = n.sqrt(sum([n.sum(n.abs(_resdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))

        # Always accept any downhill jump, advance state if bottoming out
        tol = 1 - _score/score
        if tol < 0:
            print '    Divergence'
            print '    Failed Score: %f (%5.2f%%)' % (_score, n.round(100*_score/iscore,2))
        else: 
            resdat, score, xtalk = _resdat, _score, _xtalk
            if mode == 0: srcest_bm = _srcest_bm
            elif mode == 1: srcest_ant = _srcest_ant
            else: srcest_bl = _srcest_bl

        # If we are at tolerance, add degrees of freedom and repeat
        if tol < opts.clean:
            if mode < FINAL_MODE:
                mode += 1
            else:
                if dw == opts.dw and drw == opts.drw:
                    if tier+1 < len(srctier):
                        tier += 1
                        mode = 0
                        dw,drw = 0,0
                    else:
                        print '    Tolerance'
                        break
                else:
                    mode = 0
                    dw += 1; drw += 1
        print '    New Score: %f (%5.2f%%, tol=%f)' % (score, n.round(100*score/iscore,2), tol)
    print '    Final Score: %f (%5.2f%%, tol=%f)' % (score, n.round(100*score/iscore,2), tol)

    # Pare down data and save to file
    i0 = len(times) / 6
    times = times[i0:-i0]
    n.savez('%s__times.npz' % (arg), times=n.array(times))
    n.savez('%s__afreqs.npz' % (arg), freqs=afreqs)
    for k in srcest_bm: srcest_bm[k] = srcest_bm[k][i0:-i0]
    n.savez( '%s__srcest_bm.npz' % (arg), **srcest_bm)
    __xtalk = {}
    for bl in xtalk: __xtalk[str(bl)] = xtalk[bl]
    n.savez( '%s__xtalk.npz' % (arg), **__xtalk)
    if FINAL_MODE >= 1:
        for k in srcest_ant:
            d = {}
            for i in srcest_ant[k]:
                d[str(i)] = srcest_ant[k][i][i0:-i0]
            n.savez( '%s__srcest_ant__%s.npz' % (arg,k), **d)
    if FINAL_MODE >= 2:
        for k in srcest_bl:
            d = {}
            for bl in srcest_bl[k]:
                d[str(bl)] = srcest_bl[k][bl][i0:-i0]
            n.savez( '%s__srcest_bl__%s.npz' % (arg,k), **d)

