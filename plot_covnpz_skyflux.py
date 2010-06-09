#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a
import optparse

colors = 'kbrgcmy'

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-q', '--quality', dest='quality', type='float', default=0.5,
    help='Cutoff level in quality of correlation between measured and predicted profiles for plotting sources')
opts, args = o.parse_args(sys.argv[1:])

qualities = {
    '0:42:53.96_52:07:34.8': 0.84,
   '10:01:31.41_28:48:04.0': 0.90,
   '11:14:38.91_40:37:12.7': 0.90,
   '14:11:21.08_52:07:34.8': 0.93,
   '15:04:55.31_26:01:38.9': 0.95,
   '16:28:35.62_39:32:51.3': 0.94,
    '16:51:05.63_5:00:17.4': 0.96,
   '17:20:37.50_-0:58:11.6': 0.89,
   '18:35:09.37_32:42:30.5': 0.84,
   '18:44:20.21_45:29:16.6': 0.68,
    '18:56:36.10_1:20:34.8': 0.84,
    '1:08:54.37_13:19:28.8': 0.83,
    '1:36:19.69_20:58:54.8': 0.82,
    '1:37:22.97_33:09:10.4': 0.84,
    '1:57:25.31_28:53:10.6': 0.75,
   '20:14:17.81_23:33:42.8': 0.91,
   '20:19:55.31_29:44:30.8': 0.80,
   '21:19:07.35_49:36:18.2': 0.70,
   '21:23:54.38_25:02:07.2': 0.88,
   '21:44:17.81_28:12:24.7': 0.74,
   '21:55:53.91_37:55:17.9': 0.84,
   '22:45:49.22_39:38:39.8': 0.92,
    '3:19:41.25_41:30:38.7': 0.90,
    '4:08:02.40_43:00:31.1': 0.87,
    '4:18:02.81_38:00:58.6': 0.87,
    '4:37:01.87_29:44:30.8': 0.92,
    '5:04:48.28_38:06:39.8': 0.91,
    '5:42:50.23_49:53:49.1': 0.88,
    '8:13:17.32_48:14:20.5': 0.93,
    '9:21:18.65_45:41:07.2': 0.91,
                      'Sun': 0.99,
                      'cas': 0.99,
                     'crab': 0.99,
                      'cyg': 1.00,
                      'vir': 0.98,
}
qsrcs = [s for s in qualities if qualities[s] >= .5]
#qsrcs = [s for s in qualities if qualities[s] >= .98]
#qsrcs = [s for s in qualities if qualities[s] > .92 and qualities[s] < 0.98]
#qsrcs = [s for s in qualities if qualities[s] > .90 and qualities[s] <= 0.92]
#qsrcs = [s for s in qualities if qualities[s] > .88 and qualities[s] <= 0.90]
print qsrcs


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
        if len(fwords) == 3 and not fwords[2] in qsrcs: continue
        print filename
        try: f = n.load(filename)
        except:
            print '    Load file failed'
            continue
        if fwords[1] == 'times': times = f['times']
        elif fwords[1] == 'afreqs': afreqs= f['freqs']
        elif fwords[1] == 'srcest_bm':
            for k in f.files:
                if not k in qsrcs: continue
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
    for k in srcs:
        if not srcdata.has_key(k): srcdata[k] = {}
        bmsqrt = n.sqrt(srcest_bm.get(k,0.))
        d = {}
        for i in srcest_ant.get(k,{}):
          for j in srcest_ant.get(k,{}):
            if j <= i: continue
            ai = srcest_ant[k][i]
            aj = srcest_ant[k][j]
            d[a.miriad.ij2bl(i,j)] = (bmsqrt + ai) * n.conj(bmsqrt + aj)
        for bl in srcest_bl.get(k,{}):
            d[bl] = d.get(bl,0.) + srcest_bl[k][bl]
        flag = False
        for bl in d:
            #srcdata[k][bl] = srcdata[k].get(bl,[]) + [d[bl]]
            srcdata[k][bl] = srcdata[k].get(bl,[]) + [n.average(n.abs(d[bl]), axis=1)]
            flag = True
        if flag: srctimes[k] = srctimes.get(k,[]) + [times]
for k in srcdata:
    srctimes[k] = n.concatenate(srctimes[k], axis=0)
    for bl in srcdata[k]:
        srcdata[k][bl] = n.concatenate(srcdata[k][bl], axis=0)
srcs = srcdata.keys(); srcs.sort()
if opts.cal != None:
#    srclist = []
#    for src in srcs:
#        radec = src.split('_')
#        if len(radec) == 2:
#            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
#        srclist.append(src)
#    cat = a.cal.get_catalog(opts.cal, srclist)
    aa = a.cal.get_aa(opts.cal, afreqs)
#else: cat = {}

#if 'cyg' in srcs: srcs = ['cyg'] + srcs
norm=1
totflux = {}
nosunflux = {}
usedsrc = {}
for cnt, k in enumerate(srcs):
    d,w = 0.,0.
    for bl in srcdata[k]:
        d += srcdata[k][bl]
        w += 1
    d /= w
    t = srctimes[k]
    #order = n.argsort(t)
    #d,t = d.take(order, axis=0), t.take(order)
    #I = 1
    #shape = (int(t.shape[0]/I), I)
    #ints = shape[0] * shape[1]
    #d,t = d[:ints], t[:ints]
    #d.shape,t.shape = shape + d.shape[1:], shape
    #d,t = n.average(d, axis=1), n.average(t, axis=1)
    d *= norm
    print k, d

    ## Calculate beam response
    #bm = []
    lsts = []
    for i, jd in enumerate(t):
        aa.set_jultime(jd)
        lst = aa.sidereal_time()
        lsts.append(lst)
        if not k in usedsrc.get(lst,[]):
            totflux[lst] = totflux.get(lst, 0) + d[i]
            if k != 'Sun':
                nosunflux[lst] = nosunflux.get(lst, 0) + d[i]
            usedsrc[lst] = usedsrc.get(lst, []) + [k]
    #    cat[k].compute(aa)
    #    bm.append(aa[0].bm_response(cat[k].get_crds('top'), pol=opts.pol)**2)
    #bm = n.array(bm).squeeze()
    #spec = n.sum(d*bm, axis=0)/n.sum(bm**2, axis=0)
    #if cnt == 0 and k == 'cyg':
    #    norm = cat['cyg'].jys / spec
    #    norm.shape = (1,norm.size)
    #    continue
    #ind, flx = n.polyfit(n.log10(afreqs/.150), n.log10(spec), deg=1)
    #
    #q = n.average((d-n.average(d))*(bm - n.average(bm))) / n.std(d) / n.std(bm)
    #print '%25s: FLX=%6.1f IND=%+4.2f Q=%+4.2f' % (k, 10**flx, n.round(ind,2), n.round(q,2))
    ##d /= bm
    ##_f = flx[1:-1] - .5 * (flx[2:] + flx[:-2])
    ##q = n.sum(n.abs(_f)) / n.sum(n.abs(flx[1:-1]))
    #if q < opts.quality: continue
    color = colors[cnt%len(colors)]
    #p.subplot(211)
    #p.semilogy(t, n.average(n.abs(d), axis=1), color+',', label=k)
    p.semilogy(lsts, d, color+',', label=k)
    p.ylim(.1,1e5)

    #p.subplot(212)
    #p.loglog(afreqs, spec, color+',', label=k)
    #p.loglog(afreqs, 10**n.polyval([ind,flx], n.log10(afreqs/.150)), color+'-', label=k)
    #p.xlim(afreqs[0], afreqs[-1])
    #p.ylim(10,1e5)

#p.subplot(211)
#p.legend(loc='best')
lsts = totflux.keys(); lsts.sort()
d = [totflux[lst] for lst in lsts]
p.semilogy(lsts, d, 'k^', label='tot')
d = [nosunflux[lst] for lst in lsts]
p.semilogy(lsts, d, 'kv', label='tot')
p.show()

