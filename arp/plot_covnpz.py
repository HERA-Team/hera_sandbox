#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a
import optparse

colors = 'kbrgcmy'

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

p.rcParams['legend.fontsize'] = 6

filegroups = {}
for filename in args:
    basefile = filename.split('__')[0]
    filegroups[basefile] = filegroups.get(basefile, []) + [filename]
srcdata, caldata, srctimes = {}, {}, {}
for basefile in filegroups.keys():
    for filename in filegroups[basefile]:
        f = n.load(filename)
        if filename.find('srctimes') != -1:
            for src in f.files: srctimes[src] = srctimes.get(src,[]) + [f[src]]
        elif filename.find('srcfreqs') != -1:
            afreqs = f['freqs']
        elif filename.find('srcdata') != -1:
            for src in f.files: srcdata[src] = srcdata.get(src,[]) + [f[src]]
        elif filename.find('caldata') != -1:
            for src in f.files: caldata[src] = caldata.get(src,[]) + [f[src]]
for src in srcdata:
    srctimes[src] = n.concatenate(srctimes[src], axis=0)
    srcdata[src] = n.concatenate(srcdata[src], axis=0)
    caldata[src] = n.concatenate(caldata[src], axis=0)
srckeys = srcdata.keys()
srckeys.sort()
if opts.cal != None: cat = a.cal.get_catalog(opts.cal, srckeys)
else: cat = {}
corr = []
for src in srckeys:
    srcdata[src] = n.array(srcdata[src]).real.clip(.1,n.Inf)
    caldata[src] = n.array(caldata[src])
    if True:
        c = n.average((srcdata[src] - n.average(srcdata[src])) * (caldata[src] - n.average(caldata[src]))) / n.std(srcdata[src]) / n.std(caldata[src])
    else: c = 1
    corr.append(c)
corr = n.array(corr)
order = n.argsort(corr)
if cat.has_key('cyg'):
    cat['cyg'].update_jys(afreqs)
    o = srckeys.index('cyg')
    flux = n.sum(srcdata['cyg']*caldata['cyg'], axis=0) / n.sum(caldata['cyg']**2, axis=0)
    norm = cat['cyg'].jys / flux
else: norm = 1.
print norm
for i, o in enumerate(order):
    src = srckeys[o]
    c = corr[o]
    clr = colors[i%len(colors)]
    #if c < .1: continue
    #if c < .65: continue
    if src not in ['cyg','cas','crab','vir','Sun']: continue
    p.subplot(211)
    #srcdata[src][:,1:-1] /= n.sqrt(srcdata[src][:,2:] * srcdata[src][:,:-2])
    #p.semilogy(srctimes[src], n.std(srcdata[src][:,1:-1], axis=1), ',', label=src)
    #srcdata[src][1:-1] -= .5*(srcdata[src][2:]+srcdata[src][:-2])
    #p.semilogy(srctimes[src][1:-1], n.abs(n.median(srcdata[src][1:-1], axis=1)), ',', label=src)
    p.semilogy(srctimes[src], n.median(srcdata[src], axis=1), clr+',', label=src)
    p.ylim(.1,1e5)
    p.xlim(0, 2*n.pi)
    p.subplot(212)
    flux = n.sum(srcdata[src]*caldata[src], axis=0) / n.sum(caldata[src]**2, axis=0) * norm
    srcpoly = n.polyfit(n.log10(afreqs), n.log10(flux), deg=1)
    print src, srcpoly
    if cat.has_key(src):
        cat[src].update_jys(afreqs)
        p.loglog(afreqs, cat[src].jys, clr+'-.')
    p.loglog(afreqs, 10**n.polyval(srcpoly, n.log10(afreqs)), clr+':')
    p.loglog(afreqs, flux, clr, label=src)
    p.xlim(.100, .200)
    #if False:
    #    nullsig = n.std(srcdata[src])
    #    sig = n.std(srcdata[src] - flux*caldata[src])
    #else:
    #    wgts = caldata[src]**2 / n.sum(caldata[src]**2)
    #    nullavg = n.sum(srcdata[src] * wgts)
    #    nullsig = n.sqrt(n.sum((srcdata[src] - nullavg)**2 * wgts))
    #    sig = n.sqrt(n.sum((srcdata[src] - flux*caldata[src])**2 * wgts))
    #print '%25s: flux=%7.1f (%6.1f), cor=%4.2f, gof=%4.2f' % ( src, flux, sig, c, 1-sig/nullsig)
    #p.loglog(flux*n.median(caldata[src], axis=1), n.median(n.real(srcdata[src]), axis=1), ',', label=src)
#p.loglog(10**n.arange(-2,5,.01), 10**n.arange(-2,5,.01), 'k:')
#p.xlim(1e-2,1e5)
#p.ylim(1e-2,1e5)
p.legend(loc='best')
p.show()

