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
    curtime, meas_data = None, {}
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
          if len(meas_data) > 0:
            aa.set_jultime(curtime)
            cat.compute(aa)
            srcs = [k for k in cat.keys() if cat[k].alt > 0]
            bms = n.array([aa[0].bm_response(cat[src].get_crds('top'), pol=opts.pol[0])[0,0]**2 for src in srcs])
            asrcs = n.array([cat[src].jys[0] for src in srcs])
            srcs = [srcs[i] for i in range(len(srcs)) if bms[i]*asrcs[i] > 30]
            #srcs = [srcs[i] for i in range(len(srcs)) if bms[i]*asrcs[i] > 10]
            srcs.sort()
            NSRCS = len(srcs)
            print curtime, aa.sidereal_time()
            cov = n.matrix(n.zeros((NSRCS,NSRCS), dtype=n.complex))
            data = {}
            for src in srcs:
                data[src] = {}
                for bl in meas_data:
                    i,j = a.miriad.bl2ij(bl)
                    data[src][bl] = n.conj(aa.gen_phs(cat[src], i, j, srcshape=cat[src].srcshape, resolve_src=True))
            for si, src1 in enumerate(srcs):
                for sj, src2 in enumerate(srcs[si:]):
                    sj += si
                    _sum,_wgt = 0., 0.
                    for bl in data[src1]:
                        _sum += data[src1][bl] * n.conj(data[src2][bl])
                        _wgt += n.abs(data[src2][bl])**2
                    _sum /= _wgt
                    cov[si,sj] = _sum
                    cov[sj,si] = n.conj(_sum)
            icov = cov.I
            bms = n.array([aa[0].bm_response(cat[src].get_crds('top'), pol=opts.pol[0])[0,0]**2 for src in srcs])
            asrcs = n.array([cat[src].jys[0] for src in srcs])
            # Clean loop to iteratively remove sources
            esrcs = 0.
            CLEAN_GAIN = .2
            CLEAN_TOL = 1e-4
            cycles = 0
            while True:
                __esrcs = n.zeros((NSRCS,), dtype=n.complex)
                for si, src in enumerate(srcs):
                    _sum, _wgt = 0., 0.
                    for bl in meas_data:
                        _sum += meas_data[bl] * n.conj(data[src][bl])
                        _wgt += n.abs(data[src][bl])**2
                    _sum /= _wgt
                    __esrcs[si] = _sum
                # Don't estimate fluxes for sources dominated by sidelobes
                _icov = n.array(n.diag(__esrcs) * icov)
                _esrcs = n.sum(_icov, axis=0)
                for si in range(NSRCS):
                    #print srcs[si], n.round(n.abs(_icov[:,si])), n.argmax(n.abs(_icov[:,si]))
                    if si != n.argmax(n.abs(_icov[:,si])): _esrcs[si] = 0
                _esrcs = CLEAN_GAIN * n.real(_esrcs)
                _meas_data = {}
                res1, res2 = 0., 0.
                nsamp = 0.
                for bl in meas_data:
                  _meas_data[bl] = meas_data[bl]
                  res1 += n.sum(n.abs(_meas_data[bl])**2)
                  for si, src in enumerate(srcs):
                    _meas_data[bl] -= _esrcs[si]*data[src][bl]
                  res2 += n.sum(n.abs(_meas_data[bl])**2)
                  nsamp += _meas_data[bl].size
                #print n.sqrt(res2 / nsamp)
                if 1 - res2/res1 < CLEAN_TOL: break
                meas_data = _meas_data
                esrcs += _esrcs
                cycles += 1
            print 'Cleaned in %d cycles. Residual = %f Jy RMS' % (cycles,n.sqrt(res1/nsamp))
                
            for src, asrc, bm, esrc in zip(srcs,asrcs,bms,esrcs):
                srcdata[src] = srcdata.get(src,[]) + [esrc]
                caldata[src] = caldata.get(src,[]) + [asrc*bm]
                srctimes[src] = srctimes.get(src,[]) + [aa.sidereal_time()]
                print '   ', src, cat[src].alt, asrc, asrc*bm, n.abs(esrc)
            if False:
                print '='*70
                print n.round(cov, 2)
                print '-'*70
                print n.round(icov, 2)
                print '='*70
                for si in range(NSRCS):
                  print srcs[si],
                  for sj in range(NSRCS):
                    print n.round(_esrcs[sj] * icov[si,sj]),
                  print '=', n.round(esrcs[si], 2)
                #p.clf()
                #p.imshow(n.log10(n.abs(cov)), vmax=0., vmin=-2, interpolation='nearest')
                #p.title(str(curtime))
                #p.colorbar(shrink=.5)
                #p.xticks(range(NSRCS), srcs, rotation='vertical')
                #p.yticks(range(NSRCS), srcs)
                #p.show()
          meas_data = {}
          curtime = t
        f = f.take(chans)[0]
        if f: continue
        d = (d.take(chans) / aa.passband(i,j))[0]
        bl = a.miriad.ij2bl(i,j)
        meas_data[bl] = d

for src in srcdata:
    p.subplot(211)
    p.semilogy(srctimes[src], n.real(srcdata[src]), ',', label=src)
    p.subplot(212)
    #p.plot(srctimes[src], n.angle(srcdata[src]), ',', label=src)
    p.loglog(caldata[src], n.real(srcdata[src]), ',', label=src)
p.legend()
p.show()
    #print cnt, arg, srcs
#p.show()
