#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a
import optparse

colors = 'kbrgcmy'

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True)
opts, args = o.parse_args(sys.argv[1:])

p.rcParams['legend.fontsize'] = 6

filegroups = {}
for filename in args:
    basefile = filename.split('__')[0]
    filegroups[basefile] = filegroups.get(basefile, []) + [filename]
srcdata,srctimes = {}, {}
basefiles = filegroups.keys(); basefiles.sort()
for basefile in basefiles:
    filenames = filegroups[basefile]; filenames.sort(); filenames.reverse()
    for filename in filenames:
        print filename
        f = n.load(filename)
        if filename.find('times') != -1:
            times = f['times']
        elif filename.find('afreqs') != -1: afreqs = f['freqs']
        elif filename.find('srcdata') != -1:
            src = filename.split('.')[-2].split('_')[-1]
            if not srcdata.has_key(src): srcdata[src] = {}
            for bl in f.files: srcdata[src][int(bl)] = srcdata[src].get(int(bl),[]) + [f[bl]]
            srctimes[src] = srctimes.get(src, []) + [times]
for src in srcdata:
    for bl in srcdata[src]:
        srcdata[src][bl] = n.concatenate(srcdata[src][bl], axis=0)
    srctimes[src] = n.concatenate(srctimes[src], axis=0)

srckeys = srcdata.keys(); srckeys.sort()
if opts.cal != None:
    srclist = []
    for src in srckeys:
        radec = src.split('_')
        if len(radec) == 2:
            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
        srclist.append(src)
    cat = a.cal.get_catalog(opts.cal, srclist)
    aa = a.cal.get_aa(opts.cal, afreqs)
else: cat = {}

srckeys = ['cyg'] + srckeys

norm=1
for i, src in enumerate(srckeys):
    nbls = len(srcdata[src].keys())
    bmdata = sum([srcdata[src][bl] for bl in srcdata[src]]) / nbls
    bmdata = n.array(bmdata).real.clip(.01,n.Inf)
    d,t = bmdata, srctimes[src]
    order = n.argsort(t)
    d,t = d.take(order, axis=0), t.take(order)
    I = 8
    shape = (int(t.shape[0]/I), I)
    ints = shape[0] * shape[1]
    d,t = d[:ints], t[:ints]
    d.shape,t.shape = shape + d.shape[1:], shape
    d,t = n.average(d, axis=1), n.average(t, axis=1)
    d *= norm

    # Calculate beam response
    bm = []
    for jd in t:
        aa.set_jultime(jd)
        cat[src].compute(aa)
        bm.append(aa[0].bm_response(cat[src].get_crds('top'), pol=opts.pol)**2)
    bm = n.array(bm).squeeze()
    spec = n.sum(d*bm, axis=0)/n.sum(bm**2, axis=0)
    if i == 0 and src == 'cyg':
        norm = cat['cyg'].jys / spec
        norm.shape = (1,norm.size)
        continue
    ind, flx = n.polyfit(n.log10(afreqs/.150), n.log10(spec), deg=1)
    
    q = n.average((d-n.average(d))*(bm - n.average(bm))) / n.std(d) / n.std(bm)
    print '%25s: FLX=%6.1f IND=%+4.2f Q=%+4.2f' % (src, 10**flx, n.round(ind,2), n.round(q,2))
    #d /= bm
    #_f = flx[1:-1] - .5 * (flx[2:] + flx[:-2])
    #q = n.sum(n.abs(_f)) / n.sum(n.abs(flx[1:-1]))
    #if q < .5: continue
    clr = colors[i%len(colors)]
    p.subplot(211)
    p.semilogy(t, n.average(d, axis=1), clr+',', label=src)
    p.ylim(.1,1e5)

    p.subplot(212)
    p.loglog(afreqs, spec, clr+',', label=src)
    p.loglog(afreqs, 10**n.polyval([ind,flx], n.log10(afreqs/.150)), clr+'-', label=src)
    p.xlim(afreqs[0], afreqs[-1])
    p.ylim(10,1e5)

p.subplot(211)
p.legend(loc='best')
p.show()
#import time; time.sleep(2)

