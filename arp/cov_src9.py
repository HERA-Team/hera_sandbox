#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, math

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, chan=True, src=True,
    dec=True, ant=True, pol=True)
o.add_option('-S', '--SRC', dest='SRC', 
    help='Tier 2 sources.')
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

srclist1, cutoff1, catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
try: srclist2, cutoff2, catalogs = a.scripting.parse_srcs(opts.SRC, opts.cat)
except(AttributeError): srclist2, cutoff2 = [], None
cat1 = a.cal.get_catalog(opts.cal, srclist1, cutoff1, catalogs)
cat2 = a.cal.get_catalog(opts.cal, srclist2, cutoff2, catalogs)
cat = a.cal.get_catalog(opts.cal, srclist1+srclist2, cutoff1, catalogs)

NANT = len(aa)
CLEAN_GAIN = .2
FINAL_MODE = 1

for filename in args:
    # Gather data
    print 'Reading', filename
    msrdat, msrval, simdat = {}, {}, {}
    blwgt = {}
    for k in cat:
        simdat[k] = {}
        blwgt[k] = {}
    times = []
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
                    #srcshape=cat[k].srcshape, resolve_src=False))
                u,v,w = aa.gen_uvw(i,j,cat[k]).squeeze()
            except(a.phs.PointingError):
                simd = n.zeros_like(d)
                u = v = w = n.zeros_like(d)
            simdat[k][bl] = simdat[k].get(bl,[]) + [simd]
            blwgt[k][bl] = f * n.sqrt(u**2 + v**2)
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
    # Create residuals
    resdat = {}
    for bl in msrdat: resdat[bl] = msrdat[bl].copy()

    for iter in xrange(opts.maxiter):
        print 'Iteration %d:    DW=%d, DRW=%d, MODE=%d, TIER=%d' % (iter, dw, drw, mode, tier)
        # Derive source estimators
        if mode == 0: # BEAMFORM
            _srcest_bm = {}
            for k in simdat:
                if tier == 0 and k in cat2: continue
                d, w = 0., 0.
                for bl in simdat[k]:
                    d += resdat[bl] * n.conj(simdat[k][bl]) * blwgt[k][bl]
                    w += n.abs(simdat[k][bl])**2 * blwgt[k][bl]
                d /= n.where(w > 0, w, 1)
                padlen1 = math.ceil(d.shape[1] * .2)
                d = n.concatenate([n.fliplr(d[:,:padlen1]),d,n.fliplr(d[:,-padlen1:])], axis=1)
                d = n.fft.ifft(d, axis=1) # Delay Transform
                padlen0 = math.ceil(d.shape[0] * .2)
                d = n.concatenate([n.flipud(d[:padlen0]),d,n.flipud(d[-padlen0:])])
                d = n.fft.ifft(d, axis=0) # Delay-Rate Transform
                # Apply DDR Filter
                x1,x2 = drw, -drw
                if x2 == 0: x2 = d.shape[0]
                y1,y2 = dw, -dw
                if y2 == 0: y2 = d.shape[1]
                d[x1+1:x2] = 0
                d[:,y1+1:y2] = 0
                d = n.fft.fft(d, axis=0) # undo Delay-Rate Transform
                d = d[padlen0:-padlen0]
                d = n.fft.fft(d, axis=1) # undo Delay Transform
                d = d[:,padlen1:-padlen1]
                _srcest_bm[k] = srcest_bm.get(k, 0.) + CLEAN_GAIN * d
        elif mode == 1: # VARY PER ANT
            _srcest_ant = {}
            for k in simdat:
                if tier == 0 and k in cat2: continue
                if not srcest_bm.has_key(k): continue
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
                    padlen1 = math.ceil(_d.shape[1] * .2)
                    _d = n.concatenate([n.fliplr(_d[:,:padlen1]),_d,n.fliplr(_d[:,-padlen1:])], axis=1)
                    _d = n.fft.ifft(d[i], axis=1) # Delay Transform
                    padlen0 = math.ceil(_d.shape[0] * .2)
                    _d = n.concatenate([n.flipud(_d[:padlen0]),_d,n.flipud(_d[-padlen0:])])
                    _d = n.fft.ifft(_d, axis=0) # Delay-Rate Transform
                    # Apply DDR Filter
                    x1,x2 = drw, -drw
                    if x2 == 0: x2 = _d.shape[0]
                    y1,y2 = dw, -dw
                    if y2 == 0: y2 = _d.shape[1]
                    _d[x1+1:x2] = 0
                    _d[:,y1+1:y2] = 0
                    _d = n.fft.fft(_d, axis=0) # undo Delay-Rate Transform
                    _d = _d[padlen0:-padlen0]
                    d[i] = n.fft.fft(_d, axis=1) # undo Delay Transform
                    _d = _d[:,padlen1:-padlen1]
                    if not srcest_ant.has_key(k): srcest_ant[k] = {}
                    if not _srcest_ant.has_key(k): _srcest_ant[k] = {}
                    _srcest_ant[k][i] = srcest_ant[k].get(i,0.) + CLEAN_GAIN * d[i]
        else: # POLISH UP PER BASELINE
            _srcest_bl = {}
            for k in simdat:
                if tier == 0 and k in cat2: continue
                for bl in msrdat:
                    sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                    d = resdat[bl] * n.conj(simdat[k][bl]) / sd
                    d = n.fft.ifft(d, axis=1) # Delay Transform
                    padlen = math.ceil(d.shape[0] * .2)
                    d = n.concatenate([n.flipud(d[:padlen]),d,n.flipud(d[-padlen:])])
                    d = n.fft.ifft(d, axis=0) # Delay-Rate Transform
                    # Apply DDR Filter
                    x1,x2 = drw, -drw
                    if x2 == 0: x2 = d.shape[0]
                    y1,y2 = dw, -dw
                    if y2 == 0: y2 = d.shape[1]
                    d[x1+1:x2] = 0
                    d[:,y1+1:y2] = 0
                    d = n.fft.fft(d, axis=0) # undo Delay-Rate Transform
                    d = d[padlen:-padlen]
                    d = n.fft.fft(d, axis=1) # undo Delay Transform
                    if not srcest_bl.has_key(k): srcest_bl[k] = {}
                    if not _srcest_bl.has_key(k): _srcest_bl[k] = {}
                    _srcest_bl[k][bl] = srcest_bl[k].get(bl,0.) + CLEAN_GAIN * d
        if True:
            for _k in simdat:
                if tier == 0 and _k in cat2: continue
                __resdat = {}
                for bl in msrdat:
                    i,j = a.miriad.bl2ij(bl)
                    __resdat[bl] = msrdat[bl].copy()
                    for k in simdat:
                        gi,gj = 0,0
                        if mode == 0 and k == _k:
                            if _srcest_bm.has_key(k):
                                gi = n.sqrt(_srcest_bm[k])
                                gj = gi.copy()
                        elif srcest_bm.has_key(k):
                            gi = n.sqrt(srcest_bm[k])
                            gj = gi.copy()
                        if mode == 1 and k == _k:
                            if _srcest_ant.get(k,{}).has_key(i): gi += _srcest_ant[k][i]
                            if _srcest_ant.get(k,{}).has_key(j): gj += _srcest_ant[k][j]
                        elif srcest_ant.has_key(k):
                            if srcest_ant[k].has_key(i): gi += srcest_ant[k][i]
                            if srcest_ant[k].has_key(j): gj += srcest_ant[k][j]
                        __resdat[bl] -= gi * n.conj(gj) * simdat[k][bl]
                        sd = n.abs(simdat[k][bl]).clip(1., n.Inf)
                        if mode == 2 and k == _k: __resdat[bl] -= _srcest_bl[k][bl] * simdat[k][bl] / sd
                        elif srcest_bl.has_key(k): __resdat[bl] -= srcest_bl.get(k,{}).get(bl,0.) * simdat[k][bl] /sd
                    __resdat[bl] *= msrval[bl] # Mask out invalid data
                __score = n.sqrt(sum([n.sum(n.abs(__resdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))
                if __score >= score: 
                    print '        Deleting', _k, __score
                    if mode == 0:
                        if srcest_bm.has_key(_k): _srcest_bm[_k] = srcest_bm[_k].copy()
                        elif _srcest_bm.has_key(_k): del(_srcest_bm[_k])
                    elif mode == 1:
                        if srcest_ant.has_key(_k): _srcest_ant[_k] = srcest_ant[_k].copy()
                        elif _srcest_ant.has_key(_k): del(_srcest_ant[_k])
                    else:
                        for bl in _srcest_bl[_k]: _srcest_bl[_k][bl] = srcest_bl.get(_k,{}).get(bl,0.)
                else:
                    print '   ', _k, __score
        _resdat = {}
        for bl in msrdat:
            i,j = a.miriad.bl2ij(bl)
            _resdat[bl] = msrdat[bl].copy()
            for k in simdat:
                if tier == 0 and k in cat2: continue
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
                if mode == 2: _resdat[bl] -= _srcest_bl[k][bl] * simdat[k][bl] / sd
                elif srcest_bl.has_key(k): _resdat[bl] -= srcest_bl[k][bl] * simdat[k][bl] /sd
            _resdat[bl] *= msrval[bl] # Mask out invalid data
        _score = n.sqrt(sum([n.sum(n.abs(_resdat[bl]**2)) for bl in msrdat]) / sum([n.sum(msrval[bl]) for bl in msrdat]))

        # Always accept any downhill jump, advance state if bottoming out
        tol = 1 - _score/score
        if tol < 0:
            print '    Divergence'
            print '    Failed Score: %f (%5.2f%%)' % (_score, n.round(100*_score/iscore,2))
        else: 
            resdat, score = _resdat, _score
            if mode == 0: srcest_bm = _srcest_bm
            elif mode == 1: srcest_ant = _srcest_ant
            else: srcest_bl = _srcest_bl
        if (mode == FINAL_MODE and tol < 10*opts.clean) or (mode != FINAL_MODE and tol < opts.clean):
            if mode < FINAL_MODE:
                mode += 1
            else:
                if dw == opts.dw and drw == opts.drw:
                    if tier == 0:
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
    n.savez('%s__times.npz' % (filename), times=n.array(times))
    n.savez('%s__afreqs.npz' % (filename), freqs=afreqs)
    n.savez( '%s__srcest_bm.npz' % (filename), **srcest_bm)
    if FINAL_MODE >= 1:
        for k in srcest_ant:
            d = {}
            for i in srcest_ant[k]:
                d[str(i)] = srcest_ant[k][i]
            n.savez( '%s__srcest_ant__%s.npz' % (filename,k), **d)
    if FINAL_MODE >= 2:
        for k in srcest_bl:
            d = {}
            for bl in srcest_bl[k]:
                d[str(bl)] = srcest_bl[k][bl]
            n.savez( '%s__srcest_bl__%s.npz' % (filename,k), **d)

