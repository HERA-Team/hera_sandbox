#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, math

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True, chan=True, 
    dec=True, ant=True, pol=True)
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

srclist, cutoff, catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
NANT = len(aa)
CLEAN_GAIN = .2
#CLEAN_GAIN = 1

for filename in args:
    # Gather data
    print 'Reading', filename
    msrdat, msrval, simdat = {}, {}, {}
    for k in cat: simdat[k] = {}
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
                    #srcshape=cat[k].srcshape, resolve_src=True))
                    srcshape=cat[k].srcshape, resolve_src=False))
            except(a.phs.PointingError): simd = n.zeros_like(d)
            simdat[k][bl] = simdat[k].get(bl,[]) + [simd]
    for bl in msrdat:
        msrdat[bl] = n.array(msrdat[bl])
        msrval[bl] = n.array(msrval[bl])
        for k in simdat.keys():
            simdat[k][bl] = n.array(simdat[k][bl])
            if n.all(simdat[k][bl] == 0): del(simdat[k])
    for bl in msrdat:
        i,j = a.miriad.bl2ij(bl)
        res = msrdat[bl].copy()
        val = msrval[bl]
        iscore = score = n.sqrt(n.sum(n.abs(msrdat[bl]**2)) / n.sum(val))
        print '(%d,%d) Initial Score: %f' % (i,j,score)
        srcest = {}
        dw,drw = 0,0
        for iter in xrange(opts.maxiter):
            print 'Iteration', iter
            print '    DW=%d, DRW=%d' % (dw, drw)
            # Derive source estimators
            _srcest = {}
            for k in simdat:
                d = res * n.conj(simdat[k][bl])
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
                _srcest[k] = srcest.get(k,0.) + d
            _res = msrdat[bl].copy()
            # Remove new model
            for k in simdat:
                print '    %s: <%f> (rms=%f)' % (k, n.median(_srcest[k]), n.sqrt(n.average(_srcest[k]**2)))
                _res -= CLEAN_GAIN * _srcest[k] * simdat[k][bl]
            _score = n.sqrt(n.sum(n.abs(_res**2)) / n.sum(val))
            print '(%d,%d) New Score: %f (%5.2f%%)' % (i,j,score, n.round(100*score/iscore,2))
            # Always accept any downhill jump, stop if bottoming out
            if _score > score:
                print 'Divergence'
                break
            elif 1 - _score/score < opts.clean:
                if dw == opts.dw and drw == opts.drw:
                    if  1-_score/score < opts.clean:
                        print 'Tolerance'
                        break
                else:
                    dw = min(dw+1, opts.dw)
                    drw = min(drw+1, opts.drw)
            res,score,srcest = _res,_score,_srcest
        print '(%d,%d) New Score: %f (%5.2f%%)' % (i,j,score, n.round(100*score/iscore,2))
        for k in simdat:
            print '    %s: <%f>' % (k, n.median(_srcest[k]))
        #import time; time.sleep(4)
    #for k in srcest:
    #    ts = srcest[k].keys()
    #    ts.sort()
    #    srctimes[k] = ts
    #    srcest[k] = n.array([srcest[k][t] for t in ts])
    #n.savez( '%s__srcdata.npz' % (filename), **srcest)
    #n.savez('%s__srctimes.npz' % (filename), **srctimes)
    #n.savez('%s__srcfreqs.npz' % (filename), freqs=afreqs)

