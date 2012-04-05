#! /usr/bin/env python
import aipy as a, numpy as n, capo as C
import pylab
import sys, optparse, os
"""
Calculates antenna based corrections to co-align redundant array data.

Aaron Parsons
27 March 2012

"""

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

o = optparse.OptionParser()
o.set_usage('red_cal.py *.uv')
o.set_description(__doc__)
o.add_option('--name',type='str',default='cal0',
    help="The name of your solution. [default=cal0]")
#o.add_option('--refant',type='str',default="0,0",
#    help="Choose antenna in zero referenced grid matrix location <row>,<col>. Don't pick bottom row or rightmost column.")
opts, args = o.parse_args(sys.argv[1:])


for filename in args:
    # XXX Currently hardcoded for PSA898
    A_ = [0,16,8,24,4,20,12,28]
    B_ = [i+1 for i in A_]
    C_ = [i+2 for i in A_]
    D_ = [i+3 for i in A_]
    ANTPOS = n.array([A_, B_, C_, D_])
    
    bls = {}
    conj = {}
    for ri in range(ANTPOS.shape[0]):
        for ci in range(ANTPOS.shape[1]):
            for rj in range(ANTPOS.shape[0]):
                for cj in range(ci,ANTPOS.shape[1]):
                    if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                    sep = '%d,%d' % (rj-ri, cj-ci)
                    bls[sep] = bls.get(sep,[]) + [(ANTPOS[ri,ci],ANTPOS[rj,cj])]
    for sep in bls.keys():
        if sep == '0,0' or len(bls[sep]) < 2: del(bls[sep])
    for sep in bls:
        conj[sep] = [i>j for i,j in bls[sep]]
    
    strbls = {}
    conj_bl = {}
    for sep in bls:
        strbls[sep] = []
        bl_list = []
        for (i,j),c in zip(bls[sep],conj[sep]):
            if c: i,j = j,i
            bl_list.append(ij2bl(i,j))
            strbls[sep].append('%d_%d' % (i,j))
            conj_bl[ij2bl(i,j)] = c
        bls[sep] = bl_list
        strbls[sep] = ','.join(strbls[sep])
        print sep, strbls[sep]
    
    uv = a.miriad.UV(sys.argv[-1])
    fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
    del(uv)
    
    #seps = ['0,1']
    seps = bls.keys()
    strbls = ','.join([strbls[sep] for sep in seps])
    print '-'*70
    print strbls
    print '-'*70
    
    times, d, f = C.arp.get_dict_of_uv_data([filename], strbls, 'xx', verbose=True)
    for bl in d:
        i,j = bl2ij(bl)
        if conj_bl.has_key(bl) and conj_bl[bl]: d[bl] = n.conj(d[bl])
        if i == 8 or j == 8: d[bl] = -d[bl]  # XXX remove this line once correct script is run
    
    w = {}
    for bl in f: w[bl] = n.logical_not(f[bl]).astype(n.float)
    
    NANT = ANTPOS.size
    ANTIND = dict(zip(ANTPOS.flatten(), n.arange(ANTPOS.size)))
    
    def cmp_bls(bl1,bl2):
        i1,j1 = bl2ij(bl1)
        i2,j2 = bl2ij(bl2)
        return cmp(i1,i2)
    calrow,calcol = 0,0
    calant = ANTPOS[calrow,calcol]
    cal_bl = {}
    for sep in seps:
        cal_bl[sep] = [bl for bl in bls[sep] if bl2ij(bl)[0] == ANTPOS[calrow,calcol]]
        if len(cal_bl[sep]) == 0: # If cal ant is unavailable, pick another one
            cal_bl[sep] = bls[sep][:]
            cal_bl[sep].sort(cmp_bls)
        cal_bl[sep] = cal_bl[sep][0]
    
    P = n.zeros((3,NANT), dtype=n.float)
    M = n.zeros((3,1), dtype=n.float)
    P[0,ANTIND[ANTPOS[calrow , calcol]]] = 1e6; M[0] = 0
    P[1,ANTIND[ANTPOS[calrow,calcol+1]]] = 1e6; M[0] = 0
    P[2,ANTIND[ANTPOS[calrow+1,calcol]]] = 1e6; M[0] = 0
    print P
    P_gain = n.zeros((1,NANT), dtype=n.float)
    M_gain = n.zeros((1,1), dtype=n.float)
    P_gain[0,ANTIND[ANTPOS[calrow , calcol]]] = 1e6; M_gain[0] = 0
        
    for sep in seps:
        cbl = cal_bl[sep]
        i0,j0 = bl2ij(cbl)
        if conj_bl.has_key(cbl) and conj_bl[cbl]: i0,j0 = j0,i0
        for bl in bls[sep]:
            if not d.has_key(bl):continue
            if bl == cbl: continue
            i,j = bl2ij(bl)
            if conj_bl.has_key(bl) and conj_bl[bl]: i,j = j,i
            Pline = n.zeros((1,NANT), dtype=n.float)
            Pline[0,ANTIND[j ]] +=  1; Pline[0,ANTIND[i ]] += -1
            Pline[0,ANTIND[j0]] += -1; Pline[0,ANTIND[i0]] +=  1
            Mline = n.zeros((1,1), dtype=n.float)
            P_gain_line = n.zeros((1,NANT), dtype=n.float)
            P_gain_line[0,ANTIND[j]] += 1; P_gain_line[0,ANTIND[i]] += 1
            P_gain_line[0,ANTIND[j0]] += -1; P_gain_line[0,ANTIND[i0]] += -1
            M_gain_line = n.zeros((1,1), dtype=n.float)
            g,tau,info = C.arp.redundant_bl_cal(d[cbl],w[cbl],d[bl],w[bl],fqs,use_offset=False)
            if n.isnan(g).sum()!=312: print bl, n.isnan(g).sum()
            gain = n.ma.log10(n.ma.median(n.ma.abs(n.ma.masked_invalid(g))))
            Mline[0,0] = tau
            M_gain_line[0,0] = gain
            P = n.append(P, Pline, axis=0)
            M = n.append(M, Mline, axis=0)
            P_gain = n.append(P_gain, P_gain_line, axis=0)
            M_gain = n.append(M_gain, M_gain_line, axis=0)
            print '%2d-%2d/%2d-%2d' % (i,j,i0,j0), ''.join(['v-0+^'[int(c)+2] for c in Pline.flatten()]), '%6.2f' % Mline[0,0],
            if info['dtau'] > .1: print '*', info['dtau']
            else: print
            print '%2d-%2d/%2d-%2d' % (i,j,i0,j0), ''.join(['v-0+^'[int(c)+2] for c in P_gain_line.flatten()]), '%6.3f' % M_gain_line[0,0], 'G'
            #pylab.subplot(211); pylab.plot(fqs, n.abs(g))
            #pylab.subplot(212); pylab.plot(fqs, n.angle(g))
    
    C_tau = n.linalg.lstsq(P,M)[0]
    C_tau.shape = ANTPOS.shape
    print n.around(C_tau,2)
    
    C_gain = n.linalg.lstsq(P_gain,M_gain)[0]
    C_gain.shape = ANTPOS.shape
    print n.around(10**C_gain,3)
    
    n.savez(filename+'%s.npz'%opts.name,C=C_tau,C_gain=C_gain,inpnums=ANTPOS,time=n.mean(times))
    
#import pylab; pylab.show()
