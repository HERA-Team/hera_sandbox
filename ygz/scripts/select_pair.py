__author__ = 'yunfanzhang'

import aipy as a, numpy as n
import export_beam, quick_sort

#round values to cell size
def rnd(val, cell, decimals=0):
    return n.around(val/cell,decimals=decimals) * cell

def pair_coarse(aa, src, times, dist, dist_ini=0):
    d = {}
    t = times[20]
    aa.set_jultime(t)
    src.compute(aa)
    if dist_ini == 0: dist_ini = dist
    nants = len(aa)
    for i in range(nants):
        for j in range(i+1,nants):
            uvw = aa.gen_uvw(i,j,src=src).flatten()
            if uvw[0] < 0: uvw = -uvw
            uvw_r = rnd(uvw, dist_ini)
            uv_r = (uvw_r[0],uvw_r[1])
            new_sample = ((i,j),t,(uvw[0],uvw[1]))
            d[uv_r] = d.get(uv_r,[]) + [new_sample]
    repbl = []
    for key in d.keys():
        repbl.append(d[key][0][0])
    #print repbl
   # for i,j in repbl:
   #     print i,j, aa.gen_uvw(i,j,src=src)[0][0][0],aa.gen_uvw(i,j,src=src)[1][0][0]
    d = {}
    for t in times:
        aa.set_jultime(t)
        src.compute(aa)
        for bl in repbl:
            uvw = aa.gen_uvw(*bl,src=src).flatten()
            if uvw[0] < 0: uvw = -uvw
            uvw_r = rnd(uvw, dist)
            uv_r = (uvw_r[0],uvw_r[1])
            new_sample = (bl,t,(uvw[0],uvw[1]))
            try: samples = d[uv_r]
            except(KeyError): d[uv_r] = [new_sample]
            else:
                if samples[-1][0] == bl: continue # bail if repeat entry of same baseline
                d[uv_r].append(new_sample)
    for key in d.keys(): # remove entries with no redundancy
        if len(d[key]) < 2: del d[key]
    return d


def pair_sort(pairings, freq, fbmamp, cutoff=9000.):
    sorted = []
    for key in pairings:
        L = len(pairings[key])
        for i in range(L):  # get the points pairwise
            for j in range(i+1,L):
                pt1,pt2 = pairings[key][i],pairings[key][j]
                duv = tuple(x - y for x,y in zip(pt1[2], pt2[2]))
                val = export_beam.get_overlap(freq,fbmamp,*duv)
                if n.abs(val) > cutoff:
                    sorted.append((val,(pt1[0],pt1[1]),(pt2[0],pt2[1])))
    quick_sort.quick_sort(sorted,0,len(sorted)-1)
    return sorted

def get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2):
    aa.set_jultime(t1)
    src.compute(aa)
    uvw1 = aa.gen_uvw(*bl1,src=src).flatten()
    if uvw1[0] < 0: uvw1 = -uvw1
    aa.set_jultime(t2)
    src.compute(aa)
    uvw2 = aa.gen_uvw(*bl2,src=src).flatten()
    if uvw2[0] < 0: uvw2 = -uvw2
    duv = (uvw1[0]-uvw2[0],uvw1[1]-uvw2[1])
    return export_beam.get_overlap(freq,fbmamp,*duv)

def pair_fin(pairs_sorted,dt, aa, src, freq,fbmamp,cutoff=9000.):
    clos_app = {}  #dictionary of closest approach points
    for k in n.arange(len(pairs_sorted)):
        #print pairs_sorted[k] #>> f1
        ckey = (pairs_sorted[k][1][0],pairs_sorted[k][2][0])
        count = clos_app.get(ckey,[])
        if count == []: clos_app[ckey] = (pairs_sorted[k][0],pairs_sorted[k][1][1],pairs_sorted[k][2][1])
    final = []
    for key in clos_app:
        #print key, clos_app[key]
        bl1,bl2 = key[0],key[1]
        t1,t2 = clos_app[key][1],clos_app[key][2]
        correlation = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)
        while correlation > cutoff:
            t1,t2 = t1+dt,t2+dt
            correlation = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)

        t1,t2 = clos_app[key][1],clos_app[key][2]
        correlation = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)
        while correlation > cutoff:
            t1,t2 = t1-dt,t2-dt
            correlation = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)
