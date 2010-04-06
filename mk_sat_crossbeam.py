#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True, dec=True)
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, n.array([.1375]))

cross = {}
times = {}
for filename in args:
    print 'Reading', filename
    src_id = int(filename.split('.')[-1].split('-')[-1])
    src_name = 'fm%02d' % src_id
    src = a.cal.get_catalog(opts.cal, [src_name]).values()[0]
    if not times.has_key(src_name): times[src_name] = {}
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, 'cross', opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if not times[src_name].has_key(t):
            aa.set_jultime(t)
            src.compute(aa)
            if src.alt < 0: continue
            times[src_name][t] = {'top':src.get_crds('top', ncrd=3)}
        bl = a.miriad.ij2bl(i,j)
        times[src_name][t][bl] = d[0]
        cross[bl] = None

hmaps = {}

def triplets(L):
    for i in L:
        for j in L:
            for k in L:
                if i != j and i != k and j < k: yield i,j,k
    return

bls = cross.keys()
antlist = {}
for bl in bls:
    i,j = a.miriad.bl2ij(bl)
    antlist[i] = antlist[j] = None
antlist = antlist.keys()

for src_name in times:
    print 'Gridding', src_name
    tlist = times[src_name].keys()
    for t in tlist:
        for i in bls:
            if not times[src_name][t].has_key(i):
                del(times[src_name][t])
                break
    tlist = times[src_name].keys()
    tlist.sort()
    cross = {}
    valid = 1
    for i in bls:
        cross[i] = n.array([times[src_name][t][i] for t in tlist])
        cross[i] = n.abs(cross[i])
        ai_valid = n.where(cross[i] > 0, 1, 0)
        valid = n.logical_and(valid, ai_valid)
    top = n.array([times[src_name][t]['top'] for t in tlist])
    tx,ty,tz = top.transpose()
    tx = tx.compress(valid); ty = ty.compress(valid); tz = tz.compress(valid)
    for i in bls: cross[i] = cross[i].compress(valid)
    
    # Use amp closure to infer "autos"
    autos = {}
    for i,j,k in triplets(antlist):
        if i < j: blij = a.miriad.ij2bl(i,j)
        else: blij = a.miriad.ij2bl(j,i)
        if i < k: blik = a.miriad.ij2bl(i,k)
        else: blik = a.miriad.ij2bl(k,i)
        bljk = a.miriad.ij2bl(j,k)
        autos[i] = autos.get(i,[]) + [cross[blij]*cross[blik]/cross[bljk]]
    for i in autos: autos[i] = sum(autos[i]) / len(autos[i])
            
    a0 = sum(autos.values()) / len(autos)
    a0g = a0.sum()
    # Divide out an average gain for this antenna
    for i in autos: autos[i] /= autos[i].sum() / a0g
    # Recompute average after compensating for gain of each ant
    a0 = sum(autos.values()) / len(autos)
    a0g = a0.sum()
    valid = 1
    for i in autos:
        autos[i] /= a0
        ai = autos[i] 
        ai_valid = n.logical_and(n.where(ai > .25,1,0), n.where(ai < 4.,1,0))
        valid = n.logical_and(valid, ai_valid)
        
    tx = tx.compress(valid); ty = ty.compress(valid); tz = tz.compress(valid)
    for i in autos:
        if not hmaps.has_key(i):
            hmaps[i] = a.map.Map(nside=32)
            hmaps[i].set_interpol(False)
        ai = autos[i].compress(valid)
        hmaps[i].add((tx,ty,tz), n.ones(ai.size, dtype=n.float), ai)
    for i in cross:
        if not hmaps.has_key(i):
            hmaps[i] = a.map.Map(nside=32)
            hmaps[i].set_interpol(False)
        ai = cross[i].compress(valid)
        hmaps[i].add((tx,ty,tz), n.ones(ai.size, dtype=n.float), ai)

for i in hmaps:
    print 'Writing sat_beam_%d.fits' % (i)
    hmaps[i].to_fits('sat_beam_%d.fits' % (i))

