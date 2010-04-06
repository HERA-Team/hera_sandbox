#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True, chan=True, 
    dec=True, ant=True, pol=True)
o.add_option('-d', '--dw', dest='dw', type=int, default=5,
    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
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
    srcdata = {}
    caldata = {}
    srctimes = {}
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    curtime, meas_data, meas_valid = None, {}, {}
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
          valid = 0
          for bl in meas_valid: valid |= n.any(meas_valid[bl])
          if valid:
            aa.set_jultime(curtime)
            cat.compute(aa)
            srcs = [k for k in cat.keys() if cat[k].alt > 0]
            bms = n.array([aa[0].bm_response(cat[src].get_crds('top'), pol=opts.pol[0])[0]**2 for src in srcs])
            asrcs = n.array([cat[src].jys for src in srcs])
            #srcs = [srcs[i] for i in range(len(srcs)) if n.average(bms[i]*asrcs[i]) > 10]
            srcs.sort()
            bms = n.array([aa[0].bm_response(cat[src].get_crds('top'), pol=opts.pol[0])**2 for src in srcs])
            asrcs = n.array([cat[src].jys for src in srcs])
            NSRCS = len(srcs)
            print curtime, aa.sidereal_time()

            # Compute expected covariance of source fluxes
            cov = n.zeros((NCH,NSRCS,NSRCS), dtype=n.complex)
            data = {}
            for src in srcs:
                data[src] = {}
                for bl in meas_data:
                    i,j = a.miriad.bl2ij(bl)
                    data[src][bl] = n.conj(aa.gen_phs(cat[src], i, j, srcshape=cat[src].srcshape, resolve_src=True))
                    data[src][bl] *= meas_valid[bl] # only use valid data
            for si, src1 in enumerate(srcs):
                for sj, src2 in enumerate(srcs[si:]):
                    sj += si
                    _sum,_wgt = 0., 0.
                    for bl in data[src1]:
                        if True:
                            # Apply delay filter
                            gain = n.sqrt(n.average(meas_valid[bl]**2))
                            ker = n.fft.ifft(meas_valid[bl])
                            _d = n.fft.ifft(data[src1][bl] * n.conj(data[src2][bl]))
                            _d, info = a.deconv.clean(_d, ker, tol=opts.clean, maxiter=opts.maxiter)
                            _d += info['res'] / gain
                            y1,y2 = opts.dw, -opts.dw
                            if y2 == 0: y2 = _d.shape[0]
                            _d[y1:y2] = 0
                            _d = n.fft.fft(_d)
                        else: _d = data[src1][bl] * n.conj(data[src2][bl])
                        _sum += _d
                        _wgt += n.abs(data[src2][bl])**2 # downweight resolved data
                    _sum /= _wgt.clip(1,n.Inf)
                    cov[:,si,sj] = _sum
                    cov[:,sj,si] = n.conj(_sum)
            icov = n.zeros((NCH,NSRCS,NSRCS), dtype=n.complex)
            for ch in range(NCH):
                try: icov[ch] = n.linalg.inv(cov[ch])
                except(n.linalg.linalg.LinAlgError): pass

            # Optionally begin loop by subtracting best guess
            if True:
                esrcs = n.transpose(bms.squeeze() * asrcs)
                for bl in meas_data:
                  for si, src in enumerate(srcs):
                    meas_data[bl] -= esrcs[:,si]*data[src][bl]
            else: esrcs = 0.
            # Clean loop to iteratively remove sources
            CLEAN_GAIN = .2
            #CLEAN_TOL = 1e-4
            CLEAN_TOL = 1e-6
            cycles = 0
            mode = 0
            while mode < 2:
                #print cycles, mode
                __esrcs = n.zeros((NCH,NSRCS,1), dtype=n.complex)
                for si, src in enumerate(srcs):
                    _sum, _wgt = 0., 0.
                    for bl in meas_data:
                        if mode == 0:
                            # Apply delay filter
                            _d = n.fft.ifft(meas_data[bl] * n.conj(data[src][bl]))
                            if mode == 0:
                                gain = n.sqrt(n.average(meas_valid[bl]**2))
                                ker = n.fft.ifft(meas_valid[bl])
                                _d, info = a.deconv.clean(_d, ker, tol=opts.clean, maxiter=opts.maxiter)
                                _d += info['res'] / gain
                            y1,y2 = opts.dw, -opts.dw
                            if y2 == 0: y2 = _d.shape[0]
                            _d[y1:y2] = 0
                            _d = n.fft.fft(_d)
                        else: _d = meas_data[bl] * n.conj(data[src][bl])
                        _sum += _d
                        _wgt += n.abs(data[src][bl])**2
                    _sum /= _wgt.clip(1,n.Inf)
                    __esrcs[:,si,0] = _sum
                _icov = __esrcs * icov
                _esrcs = n.sum(_icov, axis=1)
                # Don't estimate fluxes for sources dominated by sidelobes
                if False:
                    for si in range(NSRCS):
                        #if n.any(6*n.abs(_icov[:,si,si]) < n.sqrt(n.sum(n.abs(_icov[:,:,si])**2, axis=1))):
                        if n.any(n.abs(_icov[:,si,si]) < n.max(n.abs(_icov[:,:,si]), axis=1)):
                            print 'Disabling src %s on cycle %d' % (srcs[si], cycles)
                            _esrcs[:,si] = 0
                if True:
                    print 'cycle=%d mode=%d' % (cycles,mode)
                    n.set_printoptions(precision=1)
                    print map(int, n.median(n.real(esrcs), axis=0).flatten())
                    print map(int, n.median(n.real(_esrcs), axis=0).flatten())
                    print map(int, n.median(n.real(__esrcs), axis=0).flatten())
                    print '-'*70
                _esrcs = CLEAN_GAIN * n.real(_esrcs)
                _meas_data = {}
                res1, res2 = 0., 0.
                nsamp = 0.
                # Simulate and subtract visibilities for this iteration
                for bl in meas_data:
                  _meas_data[bl] = meas_data[bl].copy()
                  res1 += n.sum(n.abs(_meas_data[bl])**2 * meas_valid[bl])
                  for si, src in enumerate(srcs):
                    _meas_data[bl] -= _esrcs[:,si]*data[src][bl]
                  res2 += n.sum(n.abs(_meas_data[bl])**2 * meas_valid[bl])
                  nsamp += n.sum(meas_valid[bl])
                if cycles == 0: first_res = res1
                meas_data = _meas_data
                esrcs += _esrcs
                cycles += 1
                if 1 - n.sqrt(res2/res1) < CLEAN_TOL or cycles > opts.maxiter:
                    mode += 1
            nsamp = nsamp.clip(1,n.Inf)
            print 'Cleaned in %d cycles. Residual = %f Jy RMS\n    (initially %f Jy, terminated at %f)' % (cycles,n.sqrt(res1/nsamp), n.sqrt(first_res/nsamp), n.sqrt(res2/nsamp))
                
            if cycles == 0:
                print '    Skipping this integration'
            else:
              for si in range(NSRCS):
                src = srcs[si]
                asrc = asrcs[si]
                bm = bms[si].squeeze()
                esrc = esrcs[:,si]
                srcdata[src] = srcdata.get(src,[]) + [esrc]
                #caldata[src] = caldata.get(src,[]) + [asrc*bm]
                caldata[src] = caldata.get(src,[]) + [bm]
                srctimes[src] = srctimes.get(src,[]) + [aa.sidereal_time()]
                print '   ', src, cat[src].alt, n.median(asrc), n.median(asrc*bm), n.median(n.abs(esrc))
          curtime, meas_data, meas_valid = t, {}, {}
        bl = a.miriad.ij2bl(i,j)
        meas_data[bl] = d.take(chans) / aa.passband(i,j)
        meas_valid[bl] = n.logical_not(f.take(chans)).astype(n.int)
    n.savez(filename+'__srctimes.npz', **srctimes)
    n.savez(filename+'__srcfreqs.npz', freqs=afreqs)
    n.savez(filename+'__srcdata.npz', **srcdata)
    n.savez(filename+'__caldata.npz', **caldata)

