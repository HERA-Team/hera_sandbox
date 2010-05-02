#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True, chan=True, 
    dec=True, ant=True, pol=True)
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
NCH = chans.size
aa.select_chans(chans)

srclist, cutoff, catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
NANT = len(aa)

srcdata = {}
caldata = {}
srctimes = {}
for filename in args:
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
            #srcs = [srcs[i] for i in range(len(srcs)) if n.average(bms[i]*asrcs[i]) > 30]
            srcs = [srcs[i] for i in range(len(srcs)) if n.average(bms[i]*asrcs[i]) > 10]
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
                        _sum += data[src1][bl] * n.conj(data[src2][bl])
                        _wgt += n.abs(data[src2][bl])**2
                    _sum /= _wgt.clip(1,n.Inf)
                    cov[:,si,sj] = _sum
                    cov[:,sj,si] = n.conj(_sum)
            icov = n.zeros((NCH,NSRCS,NSRCS), dtype=n.complex)
            for ch in range(NCH):
                try: icov[ch] = n.linalg.inv(cov[ch])
                except(n.linalg.linalg.LinAlgError): pass

            # Clean loop to iteratively remove sources
            esrcs = 0.
            CLEAN_GAIN = .2
            #CLEAN_TOL = 1e-4
            CLEAN_TOL = 1e-5
            cycles = 0
            while True:
                __esrcs = n.zeros((NCH,NSRCS,1), dtype=n.complex)
                for si, src in enumerate(srcs):
                    _sum, _wgt = 0., 0.
                    for bl in meas_data:
                        _sum += meas_data[bl] * n.conj(data[src][bl]) * meas_valid[bl]
                        _wgt += n.abs(data[src][bl])**2 * meas_valid[bl]
                    _sum /= _wgt.clip(1,n.Inf)
                    __esrcs[:,si,0] = _sum
                # Don't estimate fluxes for sources dominated by sidelobes
                _icov = __esrcs * icov
                _esrcs = n.sum(_icov, axis=1)
                for si in range(NSRCS):
                    #if n.any(1.5*n.abs(_icov[:,si,si]) < n.sqrt(n.sum(n.abs(_icov[:,:,si])**2, axis=1))):
                    if n.any(6*n.abs(_icov[:,si,si]) < n.sqrt(n.sum(n.abs(_icov[:,:,si])**2, axis=1))):
                        _esrcs[:,si] = 0
                _esrcs = CLEAN_GAIN * n.real(_esrcs)
                #print n.round(_esrcs)
                _meas_data = {}
                res1, res2 = 0., 0.
                nsamp = 0.
                for bl in meas_data:
                  _meas_data[bl] = meas_data[bl].copy()
                  res1 += n.sum(n.abs(_meas_data[bl])**2 * meas_valid[bl])
                  #print '-'*70
                  #print a.miriad.bl2ij(bl)
                  #print n.sum(n.abs(_meas_data[bl])**2)
                  for si, src in enumerate(srcs):
                    #print src, _esrcs[:,si], n.abs(data[src][bl])
                    #_res1 = n.abs(_meas_data[bl])**2
                    #print _res1
                    _meas_data[bl] -= _esrcs[:,si]*data[src][bl]
                    #_res2 = n.abs(_meas_data[bl])**2
                    #print _res2
                  res2 += n.sum(n.abs(_meas_data[bl])**2 * meas_valid[bl])
                  nsamp += n.sum(meas_valid[bl])
                #print n.sqrt(res2 / res1)
                if cycles == 0: first_res = res1
                if 1 - n.sqrt(res2/res1) < CLEAN_TOL: break
                meas_data = _meas_data
                esrcs += _esrcs
                cycles += 1
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
                caldata[src] = caldata.get(src,[]) + [asrc*bm]
                srctimes[src] = srctimes.get(src,[]) + [aa.sidereal_time()]
                print '   ', src, cat[src].alt, n.median(asrc), n.median(asrc*bm), n.median(n.abs(esrc))
          curtime, meas_data, meas_valid = t, {}, {}
        bl = a.miriad.ij2bl(i,j)
        meas_data[bl] = d.take(chans) / aa.passband(i,j)
        meas_valid[bl] = n.logical_not(f.take(chans)).astype(n.int)

for src in srcdata:
    srcdata[src] = n.array(srcdata[src])
    caldata[src] = n.array(caldata[src])
    p.subplot(211)
    p.semilogy(srctimes[src], n.median(srcdata[src], axis=1), ',', label=src)
    p.ylim(.1,1e5)
    p.xlim(0, 2*n.pi)
    p.subplot(212)
    p.loglog(n.median(caldata[src], axis=1), n.median(n.real(srcdata[src]), axis=1), ',', label=src)
p.loglog(10**n.arange(-2,5,.01), 10**n.arange(-2,5,.01), 'k:')
p.xlim(1e-2,1e5)
p.ylim(1e-2,1e5)
p.legend(loc='best')
p.show()
