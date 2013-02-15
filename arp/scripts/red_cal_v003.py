#! /usr/bin/env python
"""
Calculate antenna-based corrections to co-align redundant array data.
"""

import aipy as a, numpy as n, capo
import pylab
import sys, optparse, os

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

o = optparse.OptionParser()
o.set_usage('red_cal.py *.uv')
o.set_description(__doc__)
o.add_option('--name', type='str', default='cal',
    help="Output name of solution npz. [default=cal]")
o.add_option('--plot', action='store_true',
    help="Plot the gain and phase residuals after removing the parameters solved for in this script.")
o.add_option('--verbose', action='store_true',
    help="Print a lot of stuff.")
o.add_option('--calpos', type='str', default='1,3',
    help="X,Y position (row,col) of antenna to use as the calibration reference.  Default 1,3 (antenna 25)")
o.add_option('--calpol', type='str', default='xx',
    help="Polarization to calibrate to.  Should be xx or yy.  Default xx.")
o.add_option('--maxiter', type='int', default=10,
    help="Maximum number of iterations to run in redundant calibration.  Default 10.")
#o.add_option('--refant',type='str',default="0,0",
#    help="Choose antenna in zero referenced grid matrix location <row>,<col>. Don't pick bottom row or rightmost column.")
opts, args = o.parse_args(sys.argv[1:])

# PSA-64, JD2455903...
A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]
ANTPOS_5903 = n.array([A_, B_, C_, D_])

# PSA-128, JD2456240...
A_ = [49,41,47,19,29,28,34,51]
B_ = [10, 3,25,48,24,55,27,57]
C_ = [ 9,58, 1, 4,17,13,56,59]
D_ = [22,61,35,18, 5,32,30,23]
E_ = [20,63,42,37,40,14,54,50]
F_ = [43, 2,33, 6,52, 7,12,38]
G_ = [53,21,15,16,62,44, 0,26]
H_ = [31,45, 8,11,36,60,39,46]
ANTPOS_6240 = n.array([A_, B_, C_, D_,E_,F_,G_,H_])

ANTPOS = ANTPOS_6240

for filename in args:
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
            #valid = [0,1,2,3,16,17,18,19,8,9,10,11,24,25,26,27,4,5,6,7,20,21,22,23,12,13,14]
            #valid = [0,1,2,3,16,17,18,19,12,13,14,15,28,29,30,31]
            #if not i in valid or not j in valid: continue
            bl_list.append(ij2bl(i,j))
            strbls[sep].append('%d_%d' % (i,j))
            conj_bl[ij2bl(i,j)] = c
        bls[sep] = bl_list
        strbls[sep] = ','.join(strbls[sep])
        if opts.verbose: print sep, strbls[sep]
    
    uv = a.miriad.UV(sys.argv[-1])
    fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
    del(uv)
    
    seps = ['0,1','1,1','-1,1'] #+ ['2,1', '-2,1'] + ['0,2','1,2','-1,2']
    #seps = bls.keys()
    strbls = ','.join([strbls[sep] for sep in seps])
    pols = ['xx','yy']
    #pols = ['xx']
    NPOL = len(pols)
    if opts.verbose:
        print '-'*70
        print strbls
        print '-'*70
    
    times, d, f = capo.arp.get_dict_of_uv_data([filename], strbls, ','.join(pols), verbose=True)
    for bl in d:
        i,j = bl2ij(bl)
        for pol in d[bl]:
            if conj_bl.has_key(bl) and conj_bl[bl]: d[bl][pol] = n.conj(d[bl][pol])
            #if i == 8 or j == 8: d[bl][pol] = -d[bl][pol]  # XXX remove this line once correct script is run
    
    w = {}
    for bl in f:
        w[bl] = {}
        for pol in f[bl]:
            w[bl][pol] = n.logical_not(f[bl][pol]).astype(n.float)
    
    NANT = ANTPOS.size
    ANTIND = dict(zip(ANTPOS.flatten(), n.arange(ANTPOS.size)))
    
    def cmp_bls(bl1,bl2):
        i1,j1 = bl2ij(bl1)
        i2,j2 = bl2ij(bl2)
        if i1 != i2: return cmp(i1,i2)
        else: return cmp(j1,j2)
    #calrow,calcol = 0,0
    calrow,calcol = map(int, opts.calpos.split(','))
    calant = ANTPOS[calrow,calcol]
    calpol = opts.calpol
    cal_bl = {}
    for sep in seps:
        cal_bl[sep] = [bl for bl in bls[sep] if ANTPOS[calrow,calcol] in bl2ij(bl)]
        if len(cal_bl[sep]) == 0: # If cal ant is unavailable, pick another one
            cal_bl[sep] = bls[sep][:]
        cal_bl[sep].sort(cmp_bls)
        cal_bl[sep] = cal_bl[sep][0] # XXX picks lower number baseline that involves calant, but sometimes there are 2 and this arbitrarily picks one over another without a good reason
    
    P = { # Keeps track of antennas that are being paired together in a measurement
        'phs': n.zeros((4,NPOL,NANT), dtype=n.double),
        'amp': n.zeros((2,NPOL,NANT), dtype=n.double),
    }
    M = { # Keeps track of each measurement (delay difference or gain difference)
        'phs': n.zeros((4,1), dtype=n.double),
        'amp': n.zeros((2,1), dtype=n.double),
    }
    p0 = pols.index(calpol)
    p1 = (p0 + 1) % len(pols) # For 2-pol data, this selects other pol
    # Add absolute constraints that cannot be solved for internally,
    # such as setting the delay of the reference antenna and reference baselines to 0
    P['phs'][0,p0,ANTIND[ANTPOS[calrow , calcol]]] = 1e6; M['phs'][0] = 0
    P['phs'][1,p0,ANTIND[ANTPOS[calrow,calcol+1]]] = 1e6; M['phs'][1] = 0
    P['phs'][2,p0,ANTIND[ANTPOS[calrow+1,calcol]]] = 1e6; M['phs'][2] = 0
    P['amp'][0,p0,ANTIND[ANTPOS[calrow , calcol]]] = 1e6; M['amp'][0] = 0
    # Without cross-pol data, both pols of reference ant must be manually set to fixed values.
    P['phs'][3,p1,ANTIND[ANTPOS[calrow , calcol]]] = 1e6; M['phs'][3] = 0
    P['amp'][1,p1,ANTIND[ANTPOS[calrow , calcol]]] = 1e6; M['amp'][1] = 0

    outfile = filename+'_%s.npz' % (opts.name)
    dly0,gain0 = {}, {} # Starting point for searching for solution
    if os.path.exists(outfile):
        print '      Input starting point from', outfile
        f = open(outfile)
        npz = n.load(f)
        C_phs = npz['C_phs']
        C_amp = 10**npz['C_amp']
        for pi, pol in enumerate(pols):
            if not dly0.has_key(pi): dly0[pi],gain0[pi] = {},{}
            for i,tau,g in zip(ANTPOS.flatten(), C_phs[pi].flatten(), C_amp[pi].flatten()):
                dly0[pi][i] = tau
                gain0[pi][i] = g
        f.close()
        
    for sep in seps:
        cbl = cal_bl[sep]
        i0,j0 = bl2ij(cbl)
        if conj_bl.has_key(cbl) and conj_bl[cbl]: i0,j0 = j0,i0
        for bl in bls[sep]:
            if not d.has_key(bl):continue
            for pol in d[bl]:
                if bl == cbl and pol == calpol: continue
                p = pols.index(pol)
                i,j = bl2ij(bl)
                if conj_bl.has_key(bl) and conj_bl[bl]: i,j = j,i
                _P,_M = {}, {}
                for m in ('phs','amp'):
                    _P[m] = n.zeros((1,NPOL,NANT), dtype=n.double)
                    _M[m] = n.zeros((1,1), dtype=n.double)
                _P['phs'][0, p,ANTIND[j ]] +=  1; _P['phs'][0, p,ANTIND[i ]] += -1
                _P['phs'][0,p0,ANTIND[j0]] += -1; _P['phs'][0,p0,ANTIND[i0]] +=  1
                _P['amp'][0, p,ANTIND[j ]] +=  1; _P['amp'][0, p,ANTIND[i ]] +=  1
                _P['amp'][0,p0,ANTIND[j0]] += -1; _P['amp'][0,p0,ANTIND[i0]] += -1
                try:
                    tau0 = (dly0[p][j] - dly0[p][i]) - (dly0[p0][j0] - dly0[p0][i0])
                    print (i,j),(i0,j0), pol, dly0[p][i], dly0[p][j], dly0[p0][i0], dly0[p0][j0]
                except(KeyError): tau0 = 0
                g,tau,info = capo.arp.redundant_bl_cal(d[cbl][calpol], w[cbl][calpol], d[bl][pol], w[bl][pol],
                    fqs, use_offset=False, tau=tau0, maxiter=opts.maxiter)
                print (i,j),(i0,j0), tau, tau0
                gain = n.log10(n.median(n.abs(g)))
                _M['phs'][0,0] = tau
                _M['amp'][0,0] = gain
                for m in ('phs','amp'):
                    P[m] = n.append(P[m], _P[m], axis=0)
                    M[m] = n.append(M[m], _M[m], axis=0)
                    if opts.verbose:
                        print '%2d-%2d/%2d-%2d' % (i,j,i0,j0), ''.join(['v-0+^'[int(c)+2] for c in _P[m].flatten()]), '%6.2f' % _M[m][0,0],
                        if m == 'phs':
                            if info['dtau'] > .1: print '*', info['dtau']
                            else: print
                        else: print 'G'
                if opts.plot:
                    pylab.subplot(211); pylab.plot(fqs, n.abs(g)/10**gain)
                    pylab.subplot(212); pylab.plot(fqs, n.angle(g))
    
    C = {}
    for m in P:
        x0 = P[m].shape[0]
        P[m].shape = (x0, P[m].size/x0)
        #print P[m].shape, M[m].shape
        pinv = n.linalg.pinv(P[m]) # this succeeds where lstsq fails for some reason
        C[m] = n.dot(pinv,M[m])
        #C[m] = n.linalg.lstsq(P[m],M[m])[0]
        C[m].shape = (NPOL,) + ANTPOS.shape
    print '-' * 70
    print 'Antenna Positions:'
    print ANTPOS
    print 'Pol Order:', pols
    print 'Delays:'
    print n.around(C['phs'],2)
    print 'Gains:'
    print n.around(10**C['amp'],3)
    
    print 'Writing', filename + '_%s.npz' % opts.name
    n.savez(filename+'_%s.npz' % opts.name,
        C_phs=C['phs'], C_amp=C['amp'], antpos=ANTPOS, time=n.mean(times), pols=pols)
    
if opts.plot: pylab.show()
