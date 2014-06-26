#! /usr/bin/env python
'''
Plot the position in the UV plane of redundant visibilities at
a given time and frequency index.  A single baseline should be
provided; all baselines that are redundant with the provided
baseline will be plotted.
'''
import aipy as a, numpy as n
import capo as C
import optparse, sys, os

o = optparse.OptionParser()
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
fqs = aa.get_afreqs()
del(uv)
bls,conj = C.red.group_redundant_bls(aa.ant_layout)
RLCM = C.red.RelativeLogCalMatrix(len(aa), 1)
LCM = C.red.LogCalMatrix(len(aa), 1, len(bls)) # 1 polarization
for sep in bls:
    #if sep.endswith('0'): continue
    #if not (sep.endswith('1') or sep.endswith('2')): continue
    for cnt,bl in enumerate(bls[sep]):
        #print cnt
        LCM.add_meas_record(bl, opts.pol, sep, conj[bl])
        if not sep.endswith('1'): continue
        if cnt > 0: continue
        for bl2 in bls[sep][cnt+1:]:
            RLCM.add_meas_record(bl, opts.pol, conj[bl], bl2, opts.pol, conj[bl2])
RLCM.add_constraint(aa.ant_layout[0,0], opts.pol, 0.)
RLCM.add_constraint(aa.ant_layout[0,1], opts.pol, 0.)
RLCM.add_constraint(aa.ant_layout[1,0], opts.pol, 0.)

for filename in args:
    outfile = filename+'_logcal.npz'
    if os.path.exists(outfile):
        print 'File exists:', outfile
        continue
    times, dat, flg = C.arp.get_dict_of_uv_data([filename], 'cross', opts.pol, verbose=True)
    wgt = {}
    for bl in flg:
        wgt[bl] = {}
        for pol in flg[bl]:
            wgt[bl][pol] = n.logical_not(flg[bl][pol]).astype(n.float)
            dat[bl][pol] = n.where(wgt[bl][pol], dat[bl][pol], 0)
    # Step 1: coarse relative delay solution to remove phase wraps
    dly_meas, wgt_meas = [], []
    #wraps = {}
    for cnt,(bl1,pol1,cnj1,bl2,pol2,cnj2) in enumerate(RLCM.meas_order):
        #print cnt, len(RLCM.meas_order)
        #print a.miriad.bl2ij(bl1), a.miriad.bl2ij(bl2)
        if True:
            try:
                d1,w1 = dat[bl1][pol1], wgt[bl1][pol1]
                d2,w2 = dat[bl2][pol2], wgt[bl2][pol2]
            except(KeyError):
                dly_meas.append(0)
                wgt_meas.append(0)
                continue
            if cnj1: d1 = d1.conj()
            if cnj2: d2 = d2.conj()
            g,tau,info = C.red.redundant_bl_cal(d1, w1, d2, w2, fqs, use_offset=False, tau=0, maxiter=10)
        else:
            (i1,j1),(i2,j2) = a.miriad.bl2ij(bl1), a.miriad.bl2ij(bl2)
            if cnj1: i1,j1 = j1,i1
            if cnj2: i2,j2 = j2,i2
            tau = j2 - i2 - j1 + i1
        dly_meas.append(tau)
        wgt_meas.append(1.) # XXX need to weight this for real
        #if bl1 in bls['0,1']:
        #    #dly_meas.append(tau)
        #    #phs = n.exp(-2j*n.pi*fqs*tau)
        #    #phs.shape = (1,fqs.size)
        #    #d2 *= phs / g
        #    print a.miriad.bl2ij(bl1), a.miriad.bl2ij(bl2), tau
        #    import pylab as p
        #    for plt in xrange(d1.shape[0]):
        #        g,tau,info = C.red.redundant_bl_cal(d1[plt], w1[plt], d2[plt], w2[plt], fqs, use_offset=False, tau=0, maxiter=10)
        #        g = n.average(n.abs(g))
        #        print n.around(tau,1),
        #        phs = n.exp(-2j*n.pi*fqs*tau)
        #        d2[plt] *= phs / g
        #        p.subplot(4,4,plt+1)
        #        p.plot(d1[plt].real)
        #        p.plot(d2[plt].real)
        #        p.ylim(-300,300)
        #    print
        #    p.subplot(121); C.arp.waterfall(d1, mode='phs')
        #    p.subplot(122); C.arp.waterfall(d2, mode='phs')
    #p.show()
    dly = RLCM.invert(dly_meas, einstr='pm,m')
    for tau,(bl1,pol1,cnj1,bl2,pol2,cnj2) in zip(dly_meas, RLCM.meas_order):
        #print cnt, len(RLCM.meas_order)
        (i1,j1),(i2,j2) = a.miriad.bl2ij(bl1), a.miriad.bl2ij(bl2)
        if cnj1: i1,j1 = j1,i1
        if cnj2: i2,j2 = j2,i2
        print a.miriad.bl2ij(bl1), a.miriad.bl2ij(bl2), '->', (i1,j1), (i2,j2)
        print dly[i1][pol1], dly[j1][pol1], dly[i2][pol2], dly[j2][pol2]
        print tau, ' ?= ', dly[j2][pol2] - dly[i2][pol2] - dly[j1][pol1] + dly[i1][pol1]
        if bl1 in bls['0,1']:
            try:
                d1,w1 = dat[bl1][pol1], wgt[bl1][pol1]
                d2,w2 = dat[bl2][pol2], wgt[bl2][pol2]
            except(KeyError):
                continue
            if cnj1: d1 = d1.conj()
            if cnj2: d2 = d2.conj()
            phs1 = n.exp(-2j*n.pi*fqs*(dly[j1][pol1]-dly[i1][pol1]))
            phs2 = n.exp(-2j*n.pi*fqs*(dly[j2][pol2]-dly[i2][pol2]))
            g,tau,info = C.red.redundant_bl_cal(d1, w1, d2, w2, fqs, use_offset=False, tau=0, maxiter=10)
            g = n.average(n.abs(g))
            import pylab as p
            for plt in xrange(d1.shape[0]):
                p.subplot(4,4,plt+1)
                p.plot(d1[plt].real * phs1)
                p.plot(d2[plt].real * phs2 / g)
                p.ylim(-300,300)
            #import pylab as p
            #p.subplot(121); C.arp.waterfall(d1*phs1, mode='phs')
            #p.subplot(122); C.arp.waterfall(d2*phs2, mode='phs')
            #p.show()
    p.show()
    import sys; sys.exit(0)
    # XXX apply the delay solutions, maybe subtract off average delay
    # Step 2: per-channel solution without any xtalk
    logvis = []
    for bl,pol,sep,cnj in LCM.meas_order:
        # XXX currently doesn't deal with flagged data
        i,j = a.miriad.bl2ij(bl)
        phs = n.exp(-2j*n.pi*fqs*(dly[j][opts.pol]-dly[i][opts.pol]))
        phs.shape = (1, phs.size)
        d = dat[bl][pol] * phs
        if cnj: d = d.conj()
        logvis.append(n.log(d))
        #try: logvis[-1] -= wraps[bl][pol]
        #except(KeyError): pass
    print 'Inverting'
    #vis = n.array(vis)
    logvis = n.array(logvis)
    ant_sol, sep_sol = LCM.invert(logvis, einstr='im,mtz')
    d_npz = {}
    for i in ant_sol:
        for pol in ant_sol[i]:
            print i, dly[i][pol]
            phs = n.exp(-2j*n.pi*fqs*dly[i][pol])
            phs.shape = (1,phs.size)
            d_npz['%d,%s'%(i,pol)] = ant_sol[i][pol] * phs
            #d_npz['%d,%s'%(i,pol)] = ant_sol[i][pol] * n.exp(-2j*n.pi*fqs*dly[i][pol])
            #d_npz['%d,%s'%(i,pol)] = ant_sol[i][pol] * n.exp(2j*n.pi*fqs*dly[i][pol])
    print 'Writing', outfile
    n.savez(outfile, **d_npz)
