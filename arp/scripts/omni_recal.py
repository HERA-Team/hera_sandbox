#! /usr/bin/env python
import aipy as ap, numpy as np, capo
import sys
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
    meta,gains,vismdl,xtalk,jds,lsts,freqs = capo.omni.from_npz(filename)
    mdl_bls = vismdl[DPOL]
    good_ants = gains[POL].copy()
    del(good_ants[10]) # XXX do this to test a known good antenna with a known solution
    use_bls = {}
    for gp in reds:
        gp1 = [(i,j) for i,j in gp if good_ants.has_key(i) ^ good_ants.has_key(j)]
        mdl = [bl for bl in gp if mdl_bls.has_key(bl)]
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
        gmdl = gains[POL][j][ti].conj() * mdl.conj()
        gsum[i] = gsum.get(i,0) + np.where(f, 0, d * gmdl.conj())
        gwgt[i] = gwgt.get(i,0) + np.where(f, 0, np.abs(gmdl)**2)
        #if i == 10:
        #    dmdl = gains[POL][i][ti] * gmdl
        #    plt.subplot(211)
        #    plt.plot(np.abs(d))
        #    plt.plot(np.abs(dmdl))
        #    plt.subplot(212)
        #    plt.plot(np.angle(d))
        #    plt.plot(np.angle(dmdl))
        #    plt.show()
    bad_ants = gsols[-1][0].keys()
    #for i in gains[POL]:
    #for i in [10]:
    #    plt.semilogy(np.abs(np.average(gains[POL][i], axis=0)), label=str(i), color='k')
    for i in bad_ants:
    #for i in [10]:
        gains[POL][i] = np.array([gsum[i]/gwgt[i] for gsum,gwgt in gsols])
        plt.semilogy(np.abs(np.average(gains[POL][i], axis=0)), label=str(i))#, color='r')
    plt.legend(loc='best')
    plt.show()
    import IPython; IPython.embed()
