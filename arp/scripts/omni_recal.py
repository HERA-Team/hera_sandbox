#! /usr/bin/env python
import aipy as ap, numpy as np, capo
import sys, os
import pylab as plt

POL = 'x'
DPOL = POL + POL
CAL = 'psa6622_v003'
aa = ap.cal.get_aa(CAL, np.array([.15]))
info = capo.omni.aa_to_info(aa)
reds = info.get_reds()

cal_ants = {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    outfile = filename[:-len('.npz')] + '2.npz'
    if os.path.exists(outfile):
        print '    %s exists.  Skipping...' % outfile
        continue
    meta,gains,vismdl,xtalk,jds,lsts,freqs = capo.omni.from_npz(filename)
    mdl_bls = vismdl[DPOL]
    good_ants = gains[POL].copy()
    #del(good_ants[10]) # XXX do this to test a known good antenna with a known solution
    use_bls = {}
    for gp in reds:
        gp1 = [(i,j) for i,j in gp if good_ants.has_key(i) ^ good_ants.has_key(j)]
        mdl = [bl for bl in gp if mdl_bls.has_key(bl)] # XXX I think info should be able to do this cleaner
        if len(mdl) != 1: continue
        for i,j in gp1: use_bls[(i,j)] = mdl[0]
    # Time to work with data
    uv = ap.miriad.UV(filename[:-len('.npz')])
    times,gsols = [], []
    for (_,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or times[-1] != t:
            times.append(t); ti = len(times) - 1
            gsols.append(({},{}))
            gsum,gwgt = gsols[-1]
        if use_bls.has_key((j,i)): i,j,d = j,i,d.conj() # reorient to rendunant model
        if not use_bls.has_key((i,j)): continue
        mdl = vismdl[DPOL][use_bls[(i,j)]][ti] # get the redunant model for this int
        if good_ants.has_key(i): i,j,d,mdl = j,i,d.conj(),mdl.conj() # flip to put bad ant in front
        gmdl = gains[POL][j][ti].conj() * mdl.conj() # XXX this means mdl is backward in npz files
        gsum[i] = gsum.get(i,0) + np.where(f, 0, d * gmdl.conj())
        gwgt[i] = gwgt.get(i,0) + np.where(f, 0, np.abs(gmdl)**2)
    bad_ants = gsols[-1][0].keys()
    print bad_ants
    #for i in bad_ants:
    #    gains[POL][i] = np.array([np.where(gwgt[i] > 0, gsum[i]/gwgt[i], 0) for gsum,gwgt in gsols])
    # XXX haven't filled in xtalk at all
    capo.omni.to_npz(outfile, meta, gains, vismdl, xtalk, jds, lsts, freqs, conj=False)
