#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True)
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, n.array([.1375]))

ants = {}
times = {}
for filename in args:
    print 'Reading', filename
    src_id = int(filename.split('.')[-1].split('-')[-1])
    src_name = 'fm%02d' % src_id
    src = a.cal.get_catalog(opts.cal, [src_name]).values()[0]
    if not times.has_key(src_name): times[src_name] = {}
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, 'auto', opts.pol)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if not times[src_name].has_key(t):
            aa.set_jultime(t)
            src.compute(aa)
            if src.alt < 0: continue
            times[src_name][t] = {'top':src.get_crds('top', ncrd=3)}
        times[src_name][t][i] = d[0]
        ants[i] = None

hmaps = {}
for i in ants:
    hmaps[i] = a.map.Map(nside=32)
    hmaps[i].set_interpol(False)

for src_name in times:
    print 'Gridding', src_name
    tlist = times[src_name].keys()
    for t in tlist:
        for i in hmaps:
            if not times[src_name][t].has_key(i):
                del(times[src_name][t])
                break
    tlist = times[src_name].keys()
    tlist.sort()
    for i in ants: ants[i] = n.array([times[src_name][t][i] for t in tlist])

    for i in ants: ants[i] = n.abs(ants[i])
    a0 = sum([ants[i] for i in ants]) / len(ants)
    a0g = a0.sum()
    # Divide out an average gain for this antenna
    for i in ants: ants[i] /= ants[i].sum() / a0g
    # Recompute average after compensating for gain of each ant
    a0 = sum([ants[i] for i in ants]) / len(ants)
    a0g = a0.sum()
    valid = 1
    for i in ants:
        ants[i] /= a0
        ai = ants[i] 
        ai_valid = n.logical_and(n.where(ai > .25,1,0), n.where(ai < 4.,1,0))
        valid = n.logical_and(valid, ai_valid)
        
    top = n.array([times[src_name][t]['top'] for t in tlist])
    tx,ty,tz = top.transpose()
    tx = tx.compress(valid); ty = ty.compress(valid); tz = tz.compress(valid)
    for i in hmaps:
        ai = ants[i].compress(valid)
        hmaps[i].add((tx,ty,tz), n.ones(ai.size, dtype=n.float), ai)

for i in hmaps:
    print 'Writing sat_beam_%d.fits' % (i)
    hmaps[i].to_fits('sat_beam_%d.fits' % (i))

