#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True, dec=True)
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, n.array([.1375]))

N_BLS = 6
tlist, buf = {}, {}
for filename in args:
    print 'Reading', filename
    src_id = int(filename.split('.')[-1].split('-')[-1])
    src_name = 'fm%02d' % src_id
    if not buf.has_key(src_name):
        buf[src_name] = {}
        tlist[src_name] = []
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, 'cross', opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    curtime = None
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            if curtime != None and len(buf[src_name][curtime]) == N_BLS:
                tlist[src_name].append(curtime)
            buf[src_name][t] = {}
            curtime = t
        bl = a.miriad.ij2bl(i,j)
        d = n.abs(d[0])
        if d > 0: buf[src_name][t][bl] = d
bls = buf[src_name][tlist[src_name][-1]].keys()

data,gains = {},{}
for src_name in buf:
    data[src_name] = {}
    for bl in bls:
        data[src_name][bl] = n.array([buf[src_name][t][bl] for t in tlist[src_name]])
        gains[bl] = gains.get(bl,0) + data[src_name][bl].sum()
del(buf)

for src_name in data:
    for bl in bls: data[src_name][bl] /= gains[bl]

avg_profile = {}
for src_name in data:
    buf = n.array([data[src_name][bl] for bl in bls])
    avg_profile[src_name] = n.average(buf, axis=0)
    for bl in bls: data[src_name][bl] /= avg_profile[src_name]

pxs = {}
h = a.healpix.HealpixMap(nside=32)
for src_name in data:
    print 'Gridding', src_name
    cat = a.cal.get_catalog(opts.cal, [src_name])
    def get_px(t):
        aa.set_jultime(t)
        cat.compute(aa)
        return h.crd2px(*cat.get_crds('top', ncrd=3))
    px = n.array([get_px(t) for t in tlist[src_name]]).squeeze()
    for bl in bls:
        if not pxs.has_key(bl): pxs[bl] = {}
        for p,d in zip(px, data[src_name][bl]):
            pxs[bl][p] = pxs[bl].get(p, []) + [d]
        
for bl in bls:
    h.map *= 0
    for _p in pxs[bl]:
        h.map[_p] = n.median(pxs[bl][_p])
    print 'Writing testout_%d.fits' % bl
    h.to_fits('testout_%d.fits' % bl)
        #p.semilogy([_p]*len(pxs[bl][_p]), pxs[bl][_p], 'k.')
        #p.semilogy([_p], [n.median(pxs[bl][_p])], 'ro')
        #p.semilogy([_p], [n.average(pxs[bl][_p])], 'bo')

'''
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
'''
