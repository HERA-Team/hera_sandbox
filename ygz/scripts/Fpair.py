__author__ = 'yunfanzhang'

import aipy as a, numpy as n
import Fbeam

def rnd(val, res, decimals=0):
    return n.around(val/res,decimals=decimals) * res

def pair_coarse(aa, src, times, dist=2.):
    d = {}
    t = times[0]
    aa.set_jultime(t)
    src.compute(aa)
    nants = len(aa)
    for i in range(nants):
        for j in range(i+1,nants):
            uvw = aa.gen_uvw(i,j,src=src).flatten()
            uvw_r = rnd(uvw, dist)
            uv_r = (uvw_r[0],uvw_r[1])
            new_sample = ((i,j),t,(uvw[0],uvw[1]))
            d[uv_r] = d.get(uv_r,[]) + [new_sample]
    repbl=[]
    for key in d.keys():
        repbl.append(d[key][0][0])
    #print repbl
   # for i,j in repbl:
   #     print i,j, aa.gen_uvw(i,j,src=src)[0][0][0],aa.gen_uvw(i,j,src=src)[1][0][0]
    d = {}
    #print 'd=', d
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
            if samples[-1][0] == bl: continue # bail if repeat entry of same baseline
            d[uv_r].append(new_sample)
            #if n.sqrt(u**2+v**2)<distlim and delt>1.5*dt:
            #    antd[key].append([i,j,time])
            #    #print key
            #else:
            #    antd[u1,v1]=[[i,j,time]]
    for key in d.keys(): # remove entries with no redundancy
        if len(d[key]) < 2: del d[key]
    return d


def pair_fine(pairings, freq, fbmamp):
    d = {}
    for key in pairings:
        val = Fbeam.get_overlap(freq,fbmamp,key[0],key[1])
        d[key] = [val.flatten().flatten()]
        for entry in pairings[key]:
            #print entry[0]
            #print [entry[0],entry[1]]
            d[key] = pairings.get((key),[]) + [entry[0],entry[1]]
    return d
