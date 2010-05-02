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

for filename in args:
    # Gather data
    print 'Reading', filename
    msrdat, msrval, simdat = {}, {}, {}
    uvw = {}
    curtime = None
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            msrdat[t], msrval[t], simdat[t] = {}, {}, {}
            uvw[t] = {}
            aa.set_jultime(t)
            cat.compute(aa)
            for k in cat.keys():
                if cat[k].alt > 0:
                    simdat[t][k] = {}
                    uvw[t][k] = {}
            curtime = t
        bl = a.miriad.ij2bl(i,j)
        msrval[t][bl] = n.logical_not(f.take(chans)).astype(n.int)
        msrdat[t][bl] = d.take(chans) * msrval[t][bl] / aa.passband(i,j)
        for k in simdat[t]:
            simdat[t][k][bl] = n.conj(aa.gen_phs(cat[k], i, j, srcshape=cat[k].srcshape, resolve_src=True))
            uvw[t][k][bl] = aa.gen_uvw(i, j, cat[k])

    # Compute some other quantities
    ker,gain = {}, {}
    for t in msrval:
        ker[t], gain[t] = {}, {}
        for bl in msrval[t]:
            ker[t][bl] = n.fft.ifft(msrval[t][bl])
            gain[t][bl] = n.sqrt(n.average(msrval[t][bl]**2))
            gain[t][bl] = n.where(gain[t][bl]==0,1,gain[t][bl])
    
    # Use DDR filter + bmform to get srcflux estimators vs. time/freq
    times = simdat.keys()
    times.sort()
    srcest = {}
    resdat = {}
    curres, reswgt = 0., 0.
    # Remove previous model
    for t in msrdat:
      resdat[t] = {}
      for bl in msrdat[t]:
          resdat[t][bl] = msrdat[t][bl].copy()
          resdat[t][bl] *= msrval[t][bl]
          curres += n.sum(n.abs(resdat[t][bl])**2 * msrval[t][bl])
          reswgt += n.sum(msrval[t][bl])
    curres = n.sqrt(curres / reswgt)
    print 'Initial Residual:', curres
    NITER = 40
    for iter in xrange(NITER):
      print 'Iteration', iter
      newsrcest = {}
      for k in cat:
        _srcest, _srcwgt = {}, {}
        for t in times:
            if not simdat[t].has_key(k): continue
            for bl in simdat[t][k]:
                d = resdat[t][bl] * n.conj(simdat[t][k][bl]) 
                # Delay Transform
                _d = n.fft.ifft(d)
                # don't really need to deconv within a deconv
                if False:
                    if n.any(ker[t][bl] != 0):
                        _d, info = a.deconv.clean(_d, ker[t][bl], tol=opts.clean, maxiter=opts.maxiter)
                        _d += info['res'] / n.where(gain[t][bl]==0, 1, gain[t][bl])
        
                _srcest[bl] = _srcest.get(bl, []) + [_d]
                _srcwgt[bl] = _srcwgt.get(bl, []) + [n.abs(simdat[t][k][bl])**2]
        if len(_srcest) == 0: continue
        for bl in _srcest:
            _srcest[bl], _srcwgt[bl] = n.array(_srcest[bl]), n.array(_srcwgt[bl])
        # Beamform
        bm = 0.
        bmwgt = 0.
        t = simdat.keys()[0]
        for bl in _srcest:
            u,v,w = uvw[t][k][bl].squeeze()
            bm += _srcest[bl] * n.sqrt(u**2 + v**2)
            bmwgt += _srcwgt[bl] * n.sqrt(u**2 + v**2)
        bmwgt = n.where(bmwgt == 0, 1, bmwgt)
        #print bm.shape, bmwgt.shape
        #_srcest = sum([_srcest[bl] for bl in _srcest])
        #_srcwgt = sum([n.array(_srcwgt[bl]) for bl in _srcwgt])
        #_srcwgt = n.where(_srcwgt == 0, 1, _srcwgt)
        #_srcest /= _srcwgt
        _srcest = bm / bmwgt
        # Filter
        d = _srcest
        padlen = math.ceil(d.shape[0] * .1)
        d = n.concatenate([n.flipud(d[:padlen]), d, n.flipud(d[-padlen:])])
        # Delay-Rate Transform
        _srcest = n.fft.ifft(d, axis=0)
        #for bl in _srcest:
        #    d = _srcest[bl]
        #    padlen = math.ceil(d.shape[0] * .1)
        #    d = n.concatenate([n.flipud(d[:padlen]), d, n.flipud(d[-padlen:])])
        #    # Delay-Rate Transform
        #    _srcest[bl] = n.fft.ifft(d, axis=0)
        #    # bother deconvolving this axis?
        #if opts.drw != -1:
        if iter < NITER - 5:
            #x1,x2 = opts.drw, -opts.drw
            x1,x2 = iter, -iter
            if x2 == 0: x2 = _srcest.shape[0]
            _srcest[x1+1:x2] = 0
        #if opts.dw != -1:
        if iter < NITER - 5:
            #y1,y2 = opts.dw, -opts.dw
            y1,y2 = iter, -iter
            if y2 == 0: y2 = _srcest.shape[1]
            _srcest[:,y1+1:y2] = 0
        if False:
            noise = n.random.normal(0, scale=10., size=_srcest.size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=_srcest.size))
            noise.shape = _srcest.shape
            _srcest += noise
        _srcest = n.fft.fft(_srcest, axis=0) # undo DR Transform
        # unpad the data
        _srcest = _srcest[padlen:-padlen]
        _srcest = n.fft.fft(_srcest, axis=1) # undo D Transform
        _srcest = _srcest.real
        if not srcest.has_key(k): srcest[k] = {}
        newsrcest[k] = {}
        CLEAN_GAIN = .2
        #CLEAN_GAIN = 1
        ts = [t for t in times if simdat[t].has_key(k)]
        print '    Adding %f to %s (changing by %f)' % (CLEAN_GAIN*n.median(_srcest), k, CLEAN_GAIN*n.sqrt(n.average(_srcest**2)))
        for cnt,t in enumerate(ts):
            newsrcest[k][t] = srcest[k].get(t, 0) + CLEAN_GAIN * _srcest[cnt]
      newresdat = {}
      newres, reswgt = 0., 0.
      # Remove new model
      for t in msrdat:
        newresdat[t] = {}
        for bl in msrdat[t]:
            newresdat[t][bl] = msrdat[t][bl].copy()
            for k in cat:
                if not newsrcest.has_key(k): continue
                if not newsrcest[k].has_key(t): continue
                newresdat[t][bl] -= simdat[t][k][bl] * newsrcest[k][t]
            newresdat[t][bl] *= msrval[t][bl]
            newres += n.sum(n.abs(newresdat[t][bl])**2 * msrval[t][bl])
            reswgt += n.sum(msrval[t][bl])
      newres = n.sqrt(newres / reswgt)
      print '    New Residual:', newres
      if False:
          # Decide whether to accept this jump
          SIGJY = 10.
          jump = n.log(n.random.uniform(0,1))
          #print jump, -(newres**2-curres**2)/SIGJY**2
          print jump, -((newres-700)**2-(curres-700)**2)/SIGJY**2
          #if jump < -(newres**2-curres**2)/SIGJY**2:
          if jump < -((newres-700)**2-(curres-700)**2)/SIGJY**2:
            print '   Accepting'
            resdat,curres = newresdat,newres
            srcest = newsrcest
          else:
            print '   Rejecting'
      else:
          # Always accept any downhill jump, stop if bottoming out
          if newres > curres: break
          resdat,curres = newresdat,newres
          srcest = newsrcest
    
      
    srctimes = {}
    for k in srcest:
        ts = srcest[k].keys()
        ts.sort()
        srctimes[k] = ts
        srcest[k] = n.array([srcest[k][t] for t in ts])
    n.savez( '%s__srcdata.npz' % (filename), **srcest)
    n.savez('%s__srctimes.npz' % (filename), **srctimes)
    n.savez('%s__srcfreqs.npz' % (filename), freqs=afreqs)
    #n.savez('%s__caldata.npz' % (filename), **caldata)

