__author__ = 'yunfanzhang'

import numpy as n
import export_beam, quick_sort

#round values to cell size
def rnd(val, cell, decimals=0):
    return n.around(val/cell,decimals=decimals) * cell

#coarsely determine crossings by griding the uv plane
#Format: d[ur_rounded] = [(bl,t,(u,v)),...]
def pair_coarse(aa, src, times, dist,redundant=False):
    ant_dict,repbl = {},{}
    NU,NV = len(aa.ant_layout),len(aa.ant_layout[0])
    nants = len(aa)
    print "pair_coarse: nants, NU, NV =", nants, NU, NV
    for i in range(NU):
        for j in range(NV):
            ant_dict[aa.ant_layout[i][j]] = (i,j)  #ant_dict[random ant#]=antlayoutindex
    print ant_dict
    for i in range(nants):
        for j in range(i+1,nants):
            try: dkey = (ant_dict[i][0]-ant_dict[j][0],ant_dict[i][1]-ant_dict[j][1])
            except(KeyError): dkey = (i,j)
            else:
                if dkey[0]<0: dkey = (-dkey[0],-dkey[1])
            repbl[dkey] = (i,j)
    print "pair_coarse:", len(repbl), "representative baselines, 4432 expected"
    #print repbl
  #  d = {}
  #  t = times[10]
  #  aa.set_jultime(t)
  #  src.compute(aa)
  #  if dist_ini == 0: dist_ini = dist
  #  nants = len(aa)
  #  print "dist_ini = ", dist_ini
  #  for i in range(nants):
  #      for j in range(i+1,nants):
  ##          uvw = aa.gen_uvw(i,j,src=src).flatten()
   #         if uvw[0] < 0: uvw = -uvw
  # #         uvw_r = rnd(uvw, dist_ini)
  #          uv_r = (uvw_r[0],uvw_r[1])
  #          new_sample = ((i,j),t,(uvw[0],uvw[1]))
  #          d[uv_r] = d.get(uv_r,[]) + [new_sample]
    d = {}
    for t in times:
        aa.set_jultime(t)
        src.compute(aa)
        for key in repbl:
            bl = repbl[key]
            uvw = aa.gen_uvw(*bl,src=src).flatten()
            if uvw[0] < 0: uvw = -uvw
            uvw_r = rnd(uvw, dist)
            uv_r = (uvw_r[0],uvw_r[1])
            new_sample = (bl,t,(uvw[0],uvw[1]))
            try: samples = d[uv_r]
            except(KeyError): d[uv_r] = [new_sample]
            else:
                if samples[-1][0] == bl: continue # bail if repeat entry of same baseline
                #try: delet = samples[-2][0]
                #except(IndexError): print "nothing to worry about"
                #else:
                #    if delet == bl: continue
                d[uv_r].append(new_sample)
    for key in d.keys(): # remove entries with no redundancy
        if len(d[key]) < 2: del d[key]
    return d

#sorts the given dictionary of crossings in order of decreasing correlations
#Format: sorted = [(val,(bl1,t1),(bl2,t2),(u1,v1)),...] (u1v1 used for test plots only)
def pair_sort(pairings, freq, fbmamp, cutoff=0.):
    sorted = []
    for key in pairings:
        L = len(pairings[key])
        for i in range(L):  # get the points pairwise
            for j in range(i+1,L):
                pt1,pt2 = pairings[key][i],pairings[key][j]
                duv = tuple(x - y for x,y in zip(pt1[2], pt2[2]))
                val = export_beam.get_overlap(freq,fbmamp,*duv)
                if abs(val) > cutoff:
                    sorted.append((val,(pt1[0],pt1[1]),(pt2[0],pt2[1]),pt1[2]))
    quick_sort.quick_sort(sorted,0,len(sorted)-1)
    return sorted

#get dictionary of closest approach points, works when each two tracks only cross once (satisfied in this case)
#format: clos_app[bl1,bl2] = (val, t1, t2, (u1,v1))
def get_closest(pairs_sorted):
    clos_app = {}
    for k in n.arange(len(pairs_sorted)):
        ckey = (pairs_sorted[k][1][0],pairs_sorted[k][2][0])
        count = clos_app.get(ckey,[])
        if count == []:
            clos_app[ckey] = (pairs_sorted[k][0],pairs_sorted[k][1][1],pairs_sorted[k][2][1],pairs_sorted[k][3])
    return clos_app

#Alternative way to pair_sort + get_closest, usually faster (~n vs ~nlog(n))
#format: clos_app[bl1,bl2] = (val, t1, t2, (u1,v1))
def alter_clos(pairings, freq, fbmamp, cutoff=0.):
    clos_app = {}
    print "alter_clos: len(pairings)=", len(pairings)
    for key in pairings:
        L = len(pairings[key])
        for i in range(L):  # get the points pairwise
            for j in range(i+1,L):
                pt1,pt2 = pairings[key][i],pairings[key][j]
                if pt1[0] == pt2[0]:
                    print "alter_clos: ignore self-correlating baseline: ", pt1[0]
                    continue
                duv = tuple(x - y for x,y in zip(pt1[2], pt2[2]))
                val = export_beam.get_overlap(freq,fbmamp,*duv)
                blkey = (pt1[0],pt2[0])
                #if blkey==((92,112),(0,91)) or blkey==((0,91),(92,112)): print blkey,val, duv
                clos_app[blkey] = clos_app.get(blkey,[])+[(val,pt1[1],pt2[1],pt1[2])]
    for blkey in clos_app.keys():
        N = len(clos_app[blkey])
        max,max_val = 0,0.
        for i in range(N):
            if max_val < abs(clos_app[blkey][i][0]):
                max,max_val = i, abs(clos_app[blkey][i][0])
        clos_app[blkey] = clos_app[blkey][max]
    return clos_app

#computes correlations of baselines bl1, bl2 at times t1, t2
def get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2):
    aa.set_jultime(t1)
    src.compute(aa)
    if src.alt>0:
        uvw1 = aa.gen_uvw(*bl1,src=src).flatten()
        if uvw1[0] < 0: uvw1 = -uvw1
    else: return 0  #if src below horizon, will break out of while loop
    aa.set_jultime(t2)
    src.compute(aa)
    if src.alt>0:
        uvw2 = aa.gen_uvw(*bl2,src=src).flatten()
        if uvw2[0] < 0: uvw2 = -uvw2
        duv = (uvw1[0]-uvw2[0],uvw1[1]-uvw2[1])
    else: return 0
    #print n.sqrt(duv[0]*duv[0]+duv[1]*duv[1])
    return export_beam.get_overlap(freq,fbmamp,*duv), (uvw1,uvw2)

def get_weight(aa,bl1,bl2,uvw,multweight,noiseweight):
    weight = 1.
    ant_dict = {}
    NU,NV = len(aa.ant_layout),len(aa.ant_layout[0])
    for i in range(NU):
        for j in range(NV):
            ant_dict[aa.ant_layout[i][j]] = (i,j)  #ant_dict[random ant#]=antlayoutindex
    try: multfactor = (NU-abs(ant_dict[bl1][0]-ant_dict[bl2][0]))*(NV-abs(ant_dict[bl1][1]-ant_dict[bl2][1]))
    except(KeyError): multfactor = 1
    if multweight: weight = weight*multfactor
    noisefactor = (uvw[0]*uvw[0]+uvw[1]*uvw[1]+uvw[2]*uvw[2])**(-1.5)
    if noiseweight: weight = weight*noisefactor
    return weight

# Outputs the final array of sorted pairs of points in uv space,
# spaced in time to avoid over computing information already extracted from fringe rate filtering
# format pair_fin = [(val,(bl1,t1),(bl2,t2))...]
def pair_fin(clos_app,dt, aa, src, freq,fbmamp,multweight=False,noiseweight=False,cutoff=6000.):
    final = []
    cnt, N = 0,len(clos_app)
    for key in clos_app:
        cnt = cnt+1
        if (cnt/500)*500 == cnt:
            print 'Calculating baseline pair %d out of %d:' % (cnt,N)
        bl1,bl2 = key[0],key[1]
        t1,t2 = clos_app[key][1],clos_app[key][2]
        correlation,(uvw1,uvw2) = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)
        if correlation == 0: continue
        weight = get_weight(aa,bl1,bl2,uvw1,multweight,noiseweight)
        while correlation > cutoff:
            final.append((weight*correlation,correlation,(bl1,t1,uvw1),(bl2,t2,uvw2)))
            t1,t2 = t1+dt,t2+dt
            try: correlation,(uvw1,uvw2)  = get_corr(aa, src,freq,fbmamp, t1,t2, bl1, bl2)
            except(TypeError): correlation  = 0.
            else: weight = get_weight(aa,bl1,bl2,uvw1,multweight,noiseweight)
        t1,t2 = clos_app[key][1]-dt,clos_app[key][2]-dt
        correlation,(uvw1,uvw2)  = get_corr(aa, src, freq,fbmamp, t1,t2, bl1, bl2)
        weight = get_weight(aa,bl1,bl2,uvw1[0],uvw1[1],multweight,noiseweight)
        while correlation > cutoff:
            final.append((weight*correlation,correlation,(bl1,t1,uvw1),(bl2,t2,uvw2)))
            t1,t2 = t1-dt,t2-dt
            try: correlation,(uvw1,uvw2)  = get_corr(aa, src,freq,fbmamp, t1,t2, bl1, bl2)
            except(TypeError): correlation  = 0.
            else: weight = get_weight(aa,bl1,bl2,uvw1[0],uvw1[1],multweight,noiseweight)
    quick_sort.quick_sort(final,0,len(final)-1)
    return final

#create a test sample to plot the pairs of points
def test_sample(pairs_final,cutoff=3000.):
    pairs = []
    print len(pairs_final)
    print pairs_final[0]
    bl1,bl2 = pairs_final[0][2][0],pairs_final[0][3][0]
    for entry in pairs_final:
        if (entry[2][0],entry[3][0]) != (bl1,bl2): continue
        uvw1,uvw2 = entry[2][2],entry[3][2]
        pairs.append(((uvw1[0],uvw1[1]),(uvw2[0],uvw2[1])))
    return pairs