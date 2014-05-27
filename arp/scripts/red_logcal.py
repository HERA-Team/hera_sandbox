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
del(uv)
bls,conj = C.red.group_redundant_bls(aa.ant_layout)
LCM = C.red.LogCalMatrix(len(aa), 1, len(bls)) # 1 polarization
for sep in bls:
    #if sep.endswith('0'): continue
    #if not (sep.endswith('1') or sep.endswith('2')): continue
    for bl in bls[sep]:
        LCM.add_meas_record(bl, opts.pol, sep, conj[bl])

for filename in args:
    outfile = filename+'_logcal.npz'
    if os.path.exists(outfile):
        print 'File exists:', outfile
        continue
    times, dat, flg = C.arp.get_dict_of_uv_data([filename], 'cross', opts.pol, verbose=True)
    #vis, xtalk = [], 0
    vis = []
    for bl,pol,sep,cnj in LCM.meas_order:
        # XXX currently doesn't deal with flagged data
        d = dat[bl][pol]
        if cnj: d = d.conj()
        vis.append(d)
    print 'Inverting'
    #ant_sol, sep_sol = LCM.invert(vis-xtalk, einstr='im,mtz')
    ant_sol, sep_sol = LCM.invert(n.array(vis), einstr='im,mtz')
    d_npz = {}
    for i in ant_sol:
        for pol in ant_sol[i]:
            d_npz['%d,%s'%(i,pol)] = ant_sol[i][pol]
    print 'Writing', outfile
    n.savez(outfile, **d_npz)
