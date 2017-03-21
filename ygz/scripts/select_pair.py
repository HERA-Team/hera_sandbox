__author__ = 'yunfanzhang'

import numpy as n, aipy as a
import export_beam, quick_sort
from scipy import interpolate
import pdb

#round eues to cell size
def rnd(val, cell, decimals=0):
    return n.around(val/cell,decimals=decimals) * cell

#coarsely determine crossings by griding the uv plane
#Format: d[ur_rounded] = [(bl,t,(u,v)),...]
# def coarse_128plugin(nants):
#     if nants == 128:
#         ant_dict2 = {114:(100,100),116:(100,101),117:(100,102),118:(100,103),119:(100,104),120:(100,105)}
#         ant_dict3 = {123:(101,100),124:(101,101),125:(101,102),126:(101,103),127:(101,104)}
#         f2.write(str(ant_dict2)+'\n')
#         f2.write(str(ant_dict3)+'\n')
#     if nants == 128:
#                     try: dkey = (1,ant_dict2[i][0]-ant_dict2[j][0],ant_dict2[i][1]-ant_dict2[j][1])
#                     except(KeyError):
#                         try: dkey = (2,ant_dict3[i][0]-ant_dict3[j][0],ant_dict3[i][1]-ant_dict3[j][1])
#                         except(KeyError):
#                             #pdb.set_trace()
#                             uvw = aa.gen_uvw(i,j,src=src).flatten()
#                             if uvw[0] < 0: uvw = -uvw
#                             uvw_r = rnd(uvw, add_tol)
#                             dkey = (3,uvw_r[0],uvw_r[1])
#                         else:
#                             if dkey[1]<0 or (dkey[1]==0 and dkey[2]<0): dkey = (dkey[0],-dkey[1],-dkey[2])
#                     else:
#                         if dkey[1]<0 or (dkey[1]==0 and dkey[2]<0): dkey = (dkey[0],-dkey[1],-dkey[2])
#                 else:
#     return
def pair_coarse(aa, src, times, dist,redundant=False, add_tol=0.5, northsouth=True):
    f2 = open('./redundant_bl.out', 'w')
    f2.close()
    f2 = open('./redundant_bl.out', 'a')
    ant_dict,ant_dict2,ant_dict3,repbl = {},{},{},{}
    t_ad = times[10]
    aa.set_jultime(t_ad)
    src.compute(aa)
    NU,NV = len(aa.ant_layout),len(aa.ant_layout[0])
    nants = len(aa)
    print "pair_coarse: nants, NU, NV =", nants, NU, NV
    for i in range(NU):
        for j in range(NV):
            ant_dict[aa.ant_layout[i][j]] = (i,j)  #ant_dict[random ant#]=antlayoutindex
    f2.write(str(ant_dict)+'\n')
    for i in range(nants):
        for j in range(i+1,nants):
            try: dkey = (0,ant_dict[i][0]-ant_dict[j][0],ant_dict[i][1]-ant_dict[j][1]) #(0,ns,ew)
            except(KeyError):
                #pdb.set_trace()  #all 64 antennas should be in ant_layout
                print "Bug in rep baseline, all 64 antennas should be in ant_layout"
                break
            if dkey[1] < 0 or (dkey[1] == 0 and dkey[2] < 0):
                #dkey = (dkey[0],-dkey[1],-dkey[2])
                continue
            if northsouth: repbl[dkey] = repbl.get(dkey,[]) + [(i,j)]
            else:
                if dkey[2] != 0: repbl[dkey] = repbl.get(dkey,[]) + [(i,j)] # this version excludes northsouth baselines
    print "pair_coarse:", len(repbl), "representative baselines, 4432 expected"
    #print repbl[(0,1,1)]
    #import IPython; IPython.embed()
  #  d = {}
  #
  #
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
            #print key, repbl[key]
            if key[0] == 3 and len(repbl[key]) > 1:
                print "Found simultaneously redundant baseline:", key, repbl[key]
                f2.write("Found simultaneously redundant bls:"+str(key)+str(repbl[key]))
                continue
            else: bl = repbl[key][0]
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
    f2.close()
    return d

#sorts the given dictionary of crossings in order of decreasing correlations
#Format: sorted = [(val,(bl1,t1),(bl2,t2),(u1,v1)),...] (u1v1 used for test plots only)
def pair_sort(pairings, bm_intpl, cutoff=0.):
    sorted = []
    for key in pairings:
        L = len(pairings[key])
        for i in range(L):  # get the points pairwise
            for j in range(i+1,L):
                pt1,pt2 = pairings[key][i],pairings[key][j]
                duv = tuple(x - y for x,y in zip(pt1[2], pt2[2]))
                val = export_beam.get_overlap(bm_intpl,*duv)
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
def alter_clos(pairings, bm_intpl, cutoff=0.):
    clos_app = {}
    cnt = 0
    #print "alter_clos: len(pairings)=", len(pairings)
    for key in pairings:
        cnt = cnt+1
        #if (cnt/20)*20 == cnt:
        #    print 'alter_clos: Processing key %d out of %d:' % (cnt,len(pairings))
        L = len(pairings[key])
        for i in range(L):  # get the points pairwise
            for j in range(i+1,L):
                pt1,pt2 = pairings[key][i],pairings[key][j]
                if pt1[0] == pt2[0]:
                    #print "alter_clos: ignore self-correlating baseline: ", pt1[0]
                    continue
                duv = tuple(x - y for x,y in zip(pt1[2], pt2[2]))
                val = export_beam.get_overlap(bm_intpl,*duv)
                blkey = (pt1[0],pt2[0])
                #if blkey==((92,112),(0,91)) or blkey==((0,91),(92,112)): print blkey,val, duv
                clos_app[blkey] = clos_app.get(blkey,[])+[(val,pt1[1],pt2[1],pt1[2])]
    for blkey in clos_app.keys():
        N = len(clos_app[blkey])
        #if N > 10:
        #    print "Found simultaneously redundant baseline:", blkey
        #    del clos_app[blkey]
        #    continue
        max,max_val = 0,0.
        for i in range(N):
            if max_val < abs(clos_app[blkey][i][0]):
                max,max_val = i, abs(clos_app[blkey][i][0])
        clos_app[blkey] = clos_app[blkey][max]
    return clos_app

#computes correlations of baselines bl1, bl2 at times t1, t2
def get_corr(aa, src, bm_intpl, t1,t2, bl1, bl2):
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
    return export_beam.get_overlap(bm_intpl,*duv), (uvw1,uvw2)

def get_ovlp(aa,t1,t2,rbm2interp):
    aa.set_jultime(t1)
    ra1 = aa.radec_of(0,n.pi/2)[0]
    aa.set_jultime(t2)
    dra = aa.radec_of(0,n.pi/2)[0]-ra1
    dl = n.sin(dra)
    return rbm2interp(dl,0)

def get_weight(aa,bl1,bl2,uvw,multweight,noiseweight, ovlp=1.):
    weight = ovlp
    ant_dict = {}
    num_ant = len(aa)
    NU,NV = len(aa.ant_layout),len(aa.ant_layout[0])
    for i in range(NU):
        for j in range(NV):
            ant_dict[aa.ant_layout[i][j]] = (i,j)  #ant_dict[random ant#]=antlayoutindex
    try:
        m1 = (NU-abs(ant_dict[bl1[0]][0]-ant_dict[bl1[1]][0]))*(NV-abs(ant_dict[bl1[0]][1]-ant_dict[bl1[1]][1]))
        m2 = (NU-abs(ant_dict[bl2[0]][0]-ant_dict[bl2[1]][0]))*(NV-abs(ant_dict[bl2[0]][1]-ant_dict[bl2[1]][1]))
        multfactor = m1*m2
    except(KeyError):
        if num_ant<=64:
            print "get_weight: KeyError", bl1, bl2
            return
        else: multfactor = 1
    if multweight: weight = weight*multfactor
    noisefactor = (uvw[0]*uvw[0]+uvw[1]*uvw[1]+uvw[2]*uvw[2])**(-1.5)
    if noiseweight: weight = weight*noisefactor
    return weight

# Outputs the final array of sorted pairs of points in uv space,
# spaced in time to avoid over computing information already extracted from fringe rate filtering
# format pair_fin = [(val,(bl1,t1),(bl2,t2))...]
def pair_fin(clos_app,dt, aa, src, freq,fbmamp,multweight=True,noiseweight=True,ovlpweight=True,cutoff=9000.*0.005*0.005,puv=False):
    final = []
    cnt, N = 0,len(clos_app)
    bm_intpl = export_beam.beam_interpol(freq,fbmamp,'cubic')
    if ovlpweight: #rbm2interp: FT of sq of beam
        fbm2 = n.multiply(fbmamp,fbmamp)   #element wise square for power beam
        rbm2 = n.fft.fft2(fbm2); rbm2 = n.fft.fftshift(rbm2)
        freqlm = n.fft.fftfreq(len(freq),d=(freq[1]-freq[0])); freqlm = n.fft.fftshift(freqlm)
        print "###small imaginary components are error, nothing to worry about"
        rbm2interp = interpolate.interp2d(freqlm, freqlm, rbm2, kind='cubic')
    for key in clos_app:
        cnt = cnt+1
        if (cnt/1000)*1000 == cnt:
            print 'Calculating baseline pair %d out of %d:' % (cnt,N)
        bl1,bl2 = key[0],key[1]
        t1,t2 = clos_app[key][1],clos_app[key][2]
        correlation,(uvw1,uvw2) = get_corr(aa, src, bm_intpl, t1,t2, bl1, bl2)
        if correlation == 0: continue
        if ovlpweight: ovlp = get_ovlp(aa,t1,t2,rbm2interp)
        else: ovlp = 1.
        weight = get_weight(aa,bl1,bl2,uvw1,multweight,noiseweight,ovlp)
        #while correlation > cutoff:
        if puv: final.append((weight*correlation,correlation,(bl1,t1,uvw1),(bl2,t2,uvw2)))
        else: final.append((weight*correlation,correlation,(bl1,t1),(bl2,t2)))
        #t1,t2 = t1+dt,t2+dt
        try: correlation,(uvw1,uvw2)  = get_corr(aa, src,bm_intpl, t1,t2, bl1, bl2)
        except(TypeError): correlation  = 0.
        else:
            if ovlpweight: ovlp = get_ovlp(aa,t1,t2,rbm2interp)
            else: ovlp = 1.
            weight = get_weight(aa,bl1,bl2,uvw1,multweight,noiseweight,ovlp)

    quick_sort.quick_sort(final,0,len(final)-1)
    return final

def weight_ext(dec, aa):
    l,m = n.sin(0-aa.long), n.sin(dec-aa.lat); nn = n.sqrt(1-l*l-m*m)
    ntop = n.array([l,m,nn])
    wt = export_beam.beam_real(aa[0], ntop, pol='x', sq=False)
    #print n.array(wt).shape()
    return wt[0]
def corr_ext(t1,t2,bl1,bl2,bm_intpl,aa,DEC,RA):
    CO = 0
    for dec in DEC:
        for ra in RA:
            #src = a.fit.RadioFixedBody(0, dec, mfreq=.15)
            src = a.fit.RadioFixedBody(ra, dec, mfreq=.15)
            aa.set_jultime(t1); src.compute(aa); wt1 = weight_ext(dec,aa)
            aa.set_jultime(t2); src.compute(aa); wt2 = weight_ext(dec,aa)
            weight = wt1*wt2
            c,(uvw1,uvw2) = get_corr(aa, src, bm_intpl, t1,t2, bl1, bl2)
            CO = CO + c*weight
    return CO

#already found the top pairs, this corrects for second order effects in the time stamp
def pair_ext(fin, aa, bm_intpl, ind=1,ddec=0.05,dra=0.05):
    dt = 0.000497;
    DEC = aa.lat+n.arange(-0.6,0.6,ddec)   #in radians
    RA = n.arange(-0.5,0.5,dra)
    result = []
    print fin[ind][2], fin[ind][3]
    (bl1,t1),(bl2,t2) = fin[ind][2], fin[ind][3]
    trange = n.arange(-10,10)*dt
    T1, T2 = trange+t1, trange+t2
    Cmax, t1M, t2M = 0,t1,t2
    for tt1 in T1:
        print 'tt1=',tt1
        for tt2 in T2:
            C = corr_ext(tt1,tt2,bl1,bl2,bm_intpl,aa,DEC,RA)
            if C>Cmax:
                Cmax = C; t1M,t2M=tt1,tt2
    result = (Cmax,(bl1,t1M),(bl2,t2M),t2M-t1M)
    return result





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