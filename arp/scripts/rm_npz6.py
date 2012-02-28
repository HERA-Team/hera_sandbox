#! /usr/bin/env python
import numpy as n, aipy as a
import optparse, sys, os

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
    tnpz = n.load(filename+'__srctimes.npz')
    dnpz = n.load(filename+'__srcdata.npz')
    fnpz = n.load(args[0]+'__srcfreqs.npz')
    _afreqs = fnpz['freqs']
    fq2fq = n.array([n.argmin(n.abs(_afreqs-f)).squeeze() for f in aa.get_afreqs()])
    print fq2fq
    srcdata = {}
    for src in dnpz.files:
        for cnt, jd in enumerate(tnpz[src]):
            if not srcdata.has_key(jd): srcdata[jd] = {}
            srcdata[jd][src] = dnpz[src][cnt].take(fq2fq)

    srckeys = dnpz.files; srckeys.sort()
    srclist = []
    for src in srckeys:
        radec = src.split('_')
        if len(radec) == 2:
            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
        srclist.append(src)
    cat = a.cal.get_catalog(opts.cal, srclist)

    times = srcdata.keys(); times.sort(); times = n.array(times)
    curtime, _t = None, None
    def mfunc(uv, p, d, f):
        global curtime, _t
        crd,t,(i,j) = p
        if t != curtime:
            curtime = t
            aa.set_jultime(t)
            cat.compute(aa)
            _t = times[n.argmin(n.abs(times-t))]
            #_t = t
            print t, _t
        sd = 0
        for k in cat:
            try: amp = srcdata[_t][k]
            except(KeyError): continue
            try: phs = n.conj(aa.gen_phs(cat[k],i,j, srcshape=cat[k].srcshape, resolve_src=True))
            except(a.phs.PointingError): continue
            sd += amp * phs
        #print (i,j),t, _t, n.average(n.where(f,0,n.abs(d)**2)[240:720]),
        d -= sd * aa.passband(i,j)
        #print n.average(n.where(f,0,n.abs(d)**2)[240:720])
        return p, n.where(f, 0, d), f
        
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='RMNPZ: removed sources %s\n' % (str(srckeys)))
