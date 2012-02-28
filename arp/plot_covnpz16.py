#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a
import optparse

colors = 'kbrgcmy'

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True)
opts, args = o.parse_args(sys.argv[1:])

p.rcParams['legend.fontsize'] = 6

filegroups = {}
for cnt, filename in enumerate(args):
    basefile = filename.split('__')[0]
    filegroups[basefile] = filegroups.get(basefile, []) + [filename]
srcdata, srctimes = {}, {}
basefiles = filegroups.keys(); basefiles.sort()
for basefile in basefiles:
    filenames = filegroups[basefile]; filenames.sort(); filenames.reverse()
    srcest_bm, srcest_ant, srcest_bl = {}, {}, {}
    srcs = {}
    for filename in filenames:
        fwords = filename[:-len('.npz')].split('__')
        #if len(fwords) == 3 and not fwords[2] in qsrcs: continue
        print filename
        try: f = n.load(filename)
        except:
            print '    Load file failed'
            continue
        if fwords[1] == 'times': times = f['times']
        elif fwords[1] == 'afreqs': afreqs= f['freqs']
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
    for k in srcs:
        if not srcdata.has_key(k): srcdata[k] = {}
        d = {}
        for i in srcest_ant.get(k,{}):
          for j in srcest_ant.get(k,{}):
            if j <= i: continue
            ai = srcest_ant[k][i]
            aj = srcest_ant[k][j]
            d[a.miriad.ij2bl(i,j)] = ai * n.conj(aj)
        for bl in srcest_bl.get(k,{}):
            d[bl] = d.get(bl,0.) + srcest_bl[k][bl]
        flag = False
        for bl in d:
            srcdata[k][bl] = srcdata[k].get(bl,[]) + [d[bl]]
            flag = True
        if flag: srctimes[k] = srctimes.get(k,[]) + [times]
for k in srcdata:
    srctimes[k] = n.concatenate(srctimes[k], axis=0)
    for bl in srcdata[k]:
        srcdata[k][bl] = n.concatenate(srcdata[k][bl], axis=0)
srcs = srcdata.keys(); srcs.sort()
if opts.cal != None:
    srclist = []
    for src in srcs:
        radec = src.split('_')
        if len(radec) == 2:
            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
        srclist.append(src)
    cat = a.cal.get_catalog(opts.cal, srclist)
    aa = a.cal.get_aa(opts.cal, afreqs)
else: cat = {}

if 'cyg' in srcs: srcs = ['cyg'] + srcs
norm=1
for cnt, k in enumerate(srcs):
    if k == 'vir': pass
    d,w = 0.,0.
    for bl in srcdata[k]:
        d += srcdata[k][bl]
        w += 1
    d /= w
    t = srctimes[k]
    #order = n.argsort(t)
    #d,t = d.take(order, axis=0), t.take(order)
    I = 1
    shape = (int(t.shape[0]/I), I)
    ints = shape[0] * shape[1]
    d,t = d[:ints], t[:ints]
    d.shape,t.shape = shape + d.shape[1:], shape
    d,t = n.average(d, axis=1), n.average(t, axis=1)
    d *= norm

    # Calculate beam response
    bm = []
    lsts = []
    for jd in t:
        aa.set_jultime(jd)
        lsts.append(aa.sidereal_time())
        cat[k].compute(aa)
        bm.append(aa[0].bm_response(cat[k].get_crds('top'), pol=opts.pol)**2)
    bm = n.array(bm).squeeze()
    spec = n.sum(d*bm, axis=0)/n.sum(bm**2, axis=0)
    if cnt == 0 and k == 'cyg':
        norm = cat['cyg'].jys / spec
        norm.shape = (1,norm.size)
        continue
    ind, flx = n.polyfit(n.log10(afreqs/.150), n.log10(spec), deg=1)
    
    q = n.average((d-n.average(d))*(bm - n.average(bm))) / n.std(d) / n.std(bm)
    print '%25s: FLX=%6.1f IND=%+4.2f Q=%+4.2f' % (k, 10**flx, n.round(ind,2), n.round(q,2))
    #d /= bm
    #_f = flx[1:-1] - .5 * (flx[2:] + flx[:-2])
    #q = n.sum(n.abs(_f)) / n.sum(n.abs(flx[1:-1]))
    #if q < opts.quality: continue
    color = colors[cnt%len(colors)]
    #p.subplot(211)
    #p.semilogy(t, n.average(n.abs(d), axis=1), color+',', label=k)
    p.semilogy(lsts, n.average(n.abs(d), axis=1), color+',', label=k)
    p.ylim(.1,1e5)

    #p.subplot(212)
    #p.loglog(afreqs, spec, color+',', label=k)
    #p.loglog(afreqs, 10**n.polyval([ind,flx], n.log10(afreqs/.150)), color+'-', label=k)
    #p.xlim(afreqs[0], afreqs[-1])
    #p.ylim(10,1e5)

#p.subplot(211)
#p.legend(loc='best')
p.show()

