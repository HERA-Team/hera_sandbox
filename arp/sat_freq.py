#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, 
    cal=True, src=True, dec=True)
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
srclist, cutoff, catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)

sat_freqs = [.13705, .137225, .13725, .1372875, 
    .1373125, .137435, .13746, #00.13756, 
    .1376625, .13768, .1376875, 
    .13771, .1377175, .1377375, .1378, ]

tx_width = 2.63671875e-05
afreqs = aa.get_afreqs()
#ch_wgts = [n.sinc((afreqs - f)/tx_width)**2 for f in sat_freqs]
ch_wgts = [(1-n.abs((afreqs - f)/tx_width)).clip(0,1) for f in sat_freqs]

src_dat = {}
curtime = None
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            aa.set_jultime(t)
            cat.compute(aa)
            curtime = t
        ad = n.abs(d)
        ad = n.where(ad == 0, -n.max(ad), ad)
        for src in cat.values():
            if src.alt <= 0: continue
            name = src.src_name
            #src_dat[name] = src_dat.get(name,0) + ad * n.exp(-(n.pi/2-src.alt)**2)
            src_dat[name] = src_dat.get(name,0) + ad * src.alt**2

names = src_dat.keys(); names.sort()
for name in names:
    src_dat[name] = n.array([n.sum(src_dat[name]*wgt)/n.sum(wgt**2) \
            for wgt in ch_wgts])
    rank = n.argsort(-src_dat[name])[:4]
    print name
    for r in rank:
        chs = n.argwhere(ch_wgts[r] > .01)
        print '   ', r, sat_freqs[r], src_dat[name][r], chs.squeeze(), ch_wgts[r].take(chs).squeeze()
