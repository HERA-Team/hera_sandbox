#! /usr/bin/env python
import numpy as n, aipy as a
import optparse, sys, os, glob

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    uvi = a.miriad.UV(filename)
    uvofile = filename + 'd'
    print filename, '->', uvofile
    if os.path.exists(uvofile):
        print '    File exists: skipping'
        continue
    aa = a.cal.get_aa(opts.cal, uvi['sdf'], uvi['sfreq'], uvi['nchan'])
    srcs = {}
    srcest_bl, srcest_ant, srcest_bm = {}, {}, {}
    for npzfilename in glob.glob('cov_srcdata012/%s*.npz' % filename):
        fwords = npzfilename[:-len('.npz')].split('__')
        print npzfilename
        try: f = n.load(npzfilename)
        except(IOError): continue
        if fwords[1] == 'times': times = f['times']
        elif fwords[1] == 'afreqs': _afreqs= f['freqs']
        elif fwords[1] == 'srcest_bm':
            for k in f.files:
                srcs[k] = None
                srcest_bm[k] = f[k]
        elif fwords[1] == 'srcest_ant':
            k = fwords[2]
            srcs[k] = None
            srcest_ant[k] = {}
            for i in f.files:
                srcest_ant[k][int(i)] = f[i]
        elif fwords[1] == 'srcest_bl':
            k = fwords[2]
            srcs[k] = None
            srcest_bl[k] = {}
            for bl in f.files:
                srcest_bl[k][int(bl)] = f[bl]
    srcs = srcs.keys()
    srcdata, srctimes = {}, {}
    for k in srcs:
        if not srcdata.has_key(k): srcdata[k] = {}
        bmsqrt = n.sqrt(srcest_bm.get(k,0.))
        for bl in srcest_bl[k]:
            i,j = a.miriad.bl2ij(bl)
            ai = srcest_ant.get(k,{}).get(i,0.)
            aj = srcest_ant.get(k,{}).get(j,0.)
            #d = (bmsqrt + ai) * n.conj(bmsqrt + aj) + srcest_bl[k][bl]
            d = (bmsqrt + ai) * n.conj(bmsqrt + aj)
            srcdata[k][bl] = d
        #srctimes[k] = times
    fq2fq = n.array([n.argmin(n.abs(_afreqs-f)).squeeze() for f in aa.get_afreqs()])
    print fq2fq
    #srcdata = {}
    #for k in _srcdata
    #    for cnt, jd in enumerate(srctimes[k]):
    #        if not srcdata.has_key(jd): srcdata[jd] = {}
    #        srcdata[jd][k] = _srcdata[k][cnt].take(fq2fq)

    srckeys = srcdata.keys(); srckeys.sort()
    srclist = []
    for src in srckeys:
        radec = src.split('_')
        if len(radec) == 2:
            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
        srclist.append(src)
    cat = a.cal.get_catalog(opts.cal, srclist)

    #times = srcdata.keys(); times.sort(); times = n.array(times)
    curtime, cnt = None, None
    def mfunc(uv, p, d, f):
        global curtime, cnt
        crd,t,(i,j) = p
        if t != curtime:
            curtime = t
            aa.set_jultime(t)
            cat.compute(aa)
            cnt = n.argmin(n.abs(times-t))
            print t, cnt
        sd = 0
        for k in cat:
            try: amp = srcdata[k][bl][cnt].take(fq2fq)
            except(KeyError,IndexError): continue
            try: phs = n.conj(aa.gen_phs(cat[k],i,j, srcshape=cat[k].srcshape, resolve_src=True))
            except(a.phs.PointingError): continue
            sd += srcdata[k][bl][cnt].take(fq2fq) * phs
            sd += srcest_bl[k][bl][cnt].take(fq2fq) * phs / n.abs(phs).clip(1., n.Inf)
        d -= sd * aa.passband(i,j)
        return p, n.where(f, 0, d), f
        
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='RMNPZ: removed sources %s\n' % (str(srckeys)))
