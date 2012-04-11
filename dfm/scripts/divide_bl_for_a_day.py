#! /usr/bin/env python
"""
Divide all redundant baselines of a certain separation by a fiducial baseline. Output an npz of standard deviations
(taken on the baseline axis) of ratios, averaged over 10 minutes data files.
"""

import aipy as a, numpy as n, capo as C, pylab as P
import optparse
import sys, scipy

o = optparse.OptionParser()
o.set_description(__doc__)
o.set_usage('divide_bl_for_a_day.py [opts] *.uv')
a.scripting.add_standard_options(o,cal=True)
o.add_option('-o','--out',dest='out',type='str',default='',
    help="Prefix of output npz file. Default='', which outputs 'STDs.npz'")
o.add_option('--usecal',dest='use_cal',action='store_true',
    help='Apply calibration before division.')
o.add_option('--sep',dest='sep',default='0,1',
    help="Separation distances in <row>,<col>. Default='0,1'")
opts,args = o.parse_args(sys.argv[1:])

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

USE_CAL = opts.use_cal 

A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]
ANTPOS = n.array([A_, B_, C_, D_])

bls = {}
conj = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ri,ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
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

seps=[opts.sep]
strbls = ','.join([strbls[sep] for sep in seps])

stds = {}
stds['amp'] = []
stds['phs'] = []
stds['JDs'] = []
    
for uvfile in args:     
    uv = a.miriad.UV(uvfile)
    fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
    (uvw,t,(i,j)),d = uv.read()
    del(uv)
    if USE_CAL: 
        print 'Applying calibration information'
        aa = a.cal.get_aa(opts.cal, fqs)
        aa.set_jultime(t)
    stds['JDs'].append(t)

    plot_bls = bls[seps[0]]

    times, d, f = C.arp.get_dict_of_uv_data([uvfile], strbls, 'xx', verbose=True)
    for bl in d:
        i,j = bl2ij(bl)
        if USE_CAL:
            d[bl] = aa.phs2src(d[bl], 'z', i,j)
            d[bl] /= n.median(aa.passband(i,j))
        if conj_bl.has_key(bl) and conj_bl[bl]: d[bl] = n.conj(d[bl])
    
    w = {}
    for bl in f: w[bl] = n.logical_not(f[bl]).astype(n.float)
    
    NANT = ANTPOS.size
    ANTIND = dict(zip(ANTPOS.flatten(), n.arange(ANTPOS.size)))
    
    calrow,calcol = 0,0
    calant = ANTPOS[calrow,calcol]
    cal_bl = {}
    for sep in seps: cal_bl[sep] = [bl for bl in bls[sep] if bl2ij(bl)[0] == ANTPOS[calrow,calcol]][0]
    
    flags = n.load(uvfile+'R.npz')
    f = n.zeros_like(flags['corrflags'])
    for thing in flags.files:
        f |= flags[thing]
    
    amp = []
    phs = []
    for sep in seps:
        cbl = cal_bl[sep]
        d[cbl] = n.ma.array(d[cbl],mask=f)
        i0,j0 = bl2ij(cbl)
        if conj_bl.has_key(cbl) and conj_bl[cbl]: i0,j0 = j0,i0
        for bl in bls[sep]:
            d[bl] = n.ma.array(d[bl],mask=f)
            if bl == cbl or not bl in plot_bls: continue
            i,j = bl2ij(bl)
            if conj_bl.has_key(bl) and conj_bl[bl]: i,j = j,i
            g = n.ma.mean(d[bl]/d[cbl],axis=0) 
            amp.append(n.abs(g))
            phs.append(n.ma.arctan2(g.real,g.imag))
    s_amp = n.where(amp == 0, n.NaN, n.ma.std(amp,axis=0)) 
    s_phs = n.where(amp == 0, n.NaN, n.ma.std(phs,axis=0)) 
    
    print "AMP STD:",s_amp[500]
    print "PHS STD:",s_phs[500]
    stds['amp'].append(s_amp)
    stds['phs'].append(s_phs)

outfile = opts.out+'STDs.npz'
print 'Saving to %s'%outfile
n.savez(outfile,JDs=stds['JDs'],amp=stds['amp'],phs=stds['phs'])
