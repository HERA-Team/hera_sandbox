#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True)
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))
srclist, cutoff, catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
NANT = len(aa)

for cnt, arg in enumerate(args):
    aa.set_jultime(float(arg))
    cat.compute(aa)
    srcs = [k for k in cat.keys() if cat[k].alt > 0]
    srcs.sort()
    NSRCS = len(srcs)
    cov = n.matrix(n.zeros((NSRCS,NSRCS), dtype=n.complex))
    data = {}
    for src in srcs:
        data[src] = {}
        for i in range(NANT):
            for j in range(i,NANT):
                bl = a.miriad.ij2bl(i,j)
                try: data[src][bl] = aa.gen_phs(cat[src], i, j)
                except(a.phs.PointingError): data[src][bl] = 0.
    for si, src1 in enumerate(srcs):
        for sj, src2 in enumerate(srcs[si:]):
            sj += si
            _sum,_wgt = 0., 0.
            for bl in data[src1]:
                _sum += data[src1][bl] * n.conj(data[src2][bl])
                _wgt += 1.
            _sum /= _wgt
            cov[si,sj] = _sum
            cov[sj,si] = n.conj(_sum)

    icov = cov.I
    sim_data = {}
    for src in srcs:
        for bl in data[src]:
            sim_data[bl] = sim_data.get(bl, 0.) + cat[src].jys[0]*data[src][bl]
    _esrcs = n.zeros((NSRCS,), dtype=n.complex)
    asrcs = n.array([cat[src].jys[0] for src in srcs])
    for si, src in enumerate(srcs):
        _sum, _wgt = 0., 0.
        for bl in sim_data:
            _sum += sim_data[bl] * n.conj(data[src][bl])
            _wgt += 1.
        _sum /= _wgt
        _esrcs[si] = _sum
    esrcs = n.abs(n.array(_esrcs * icov).flatten())

    print _esrcs / asrcs - 1
    print esrcs / asrcs - 1
    
    p.subplot(1, len(args), cnt+1)
    p.imshow(n.log10(n.abs(cov)), vmax=0., vmin=-2, interpolation='nearest')
    p.title(arg)
    p.colorbar(shrink=.5)
    p.xticks(range(NSRCS), srcs, rotation='vertical')
    p.yticks(range(NSRCS), srcs)
    print cnt, arg, srcs
p.show()
