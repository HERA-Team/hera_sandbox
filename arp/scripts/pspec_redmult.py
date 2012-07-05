#! /usr/bin/env python
import aipy as a, numpy as n
import capo as C
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, dec=True)
opts,args = o.parse_args(sys.argv[1:])

B = 0.008
#B = 0.04 / 3.
NTAPS = 3
WINDOW = 'blackman-harris'

# XXX Currently hardcoded for PSA898
A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]
ANTPOS = n.array([A_, B_, C_, D_])

def bl_index(bl):
    i,j = a.miriad.bl2ij(bl)
    return i * 32 + j

# Get a dict of all separations and the bls that contribute
bl2sep = {}
sep2bl = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
                if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                sep = a.miriad.ij2bl(rj-ri, cj-ci)
                i,j = ANTPOS[ri,ci], ANTPOS[rj,cj]
                bl = a.miriad.ij2bl(i,j)
                if i > j: i,j,sep = j,i,-sep
                bl2sep[bl] = sep
                sep = n.abs(sep)
                sep2bl[sep] = sep2bl.get(sep,[]) + [bl]

uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

#cen_fqs = n.arange(.115,.190,.005)
cen_fqs = n.array([.162])
kwargs = {'cen_fqs':cen_fqs,'B':B, 'ntaps':NTAPS, 'window':WINDOW, 'bm_fqs':freqs.clip(.120,.190)}
#window = a.dsp.gen_window(freqs.size, window=WINDOW)

for filename in args:
    outfile = filename + '.pspec'
    print filename,'->',outfile
    if os.path.exists(outfile):
        print '    File exists.  Skipping...'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi, override={'nchan':32}) # XXX
    uvo._wrhd('history', uvi['history'] + 'PSPEC: ' + ' '.join(sys.argv) + '\n')
    uvo.add_var('k3pk_fq', 'r')
    uvo.add_var('k3pk_wgt', 'r')
    
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
    uvi.select('decimate', opts.decimate, opts.decphs)

    Tlist,Wlist,curtime = {},{},None
    klist = {}
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        # XXX Need to deal with polarization 
        if t != curtime:
            for sep,bls in sep2bl.items():
                if len(bls) < 2: continue
                for cnt,bl0 in enumerate(bls):
                    if not Tlist.has_key(bl0) or len(Tlist[bl0]) == 0: continue
                    #i0,j0 = a.miriad.bl2ij(bl0)
                    #if 16 in [i0,j0]: print a.miriad.bl2ij(bl0), 'check'
                    #print sep, [bl for bl in bls if len(Tlist[bl]) > 0]
                    for bl1 in bls[cnt:]: # this includes "auto-pspecs"
                        if len(Tlist[bl1]) == 0: continue
                        for fq in Tlist[bl0]:
                            _Tlist = [Tlist[bl0][fq], Tlist[bl1][fq]]
                            _Wlist = [Wlist[bl0][fq], Wlist[bl1][fq]]
                            D, wgt = C.pspec.k3pk_from_Trms(_Tlist, _Wlist, k=klist[bl0][fq][0], B=B, fq=fq)
                            D,wgt = D[0], wgt[0]
                            D /= wgt # Write everything to file in mk^2 units
                            uvo['k3pk_fq'] = fq
                            uvo['k3pk_wgt'] = wgt
                            #uvo.copyvr(uvi)
                            print a.miriad.bl2ij(bl0), a.miriad.bl2ij(bl1), t, wgt
                            uvo.write((crd,curtime,(bl_index(bl0),bl_index(bl1))), D, n.zeros(D.shape, dtype=n.int))
            # Clear the current pspec data and start a new integration
            Tlist,Wlist = {},{}
            for bl in bl2sep: Tlist[bl],Wlist[bl] = {},{}
            curtime = t

        bl = a.miriad.ij2bl(i,j)
        sep = bl2sep[bl]
        if sep < 0: d,sep = n.conj(d),-sep
        w = n.logical_not(f).astype(n.float)
        Tres,ks = C.pspec.Trms_vs_fq(freqs, d, **kwargs)
        W,ks = C.pspec.Trms_vs_fq(freqs, w, **kwargs)

        for fq in Tres:
            gain = n.abs(W[fq][0])
            print fq, (i,j), sep, t, gain
            T_cl, info = a.deconv.clean(Tres[fq], W[fq], tol=1e-9, maxiter=100, stop_if_div=False, verbose=False)
            #print '   ', int(n.around(fq*1e3)), info['term']
            T = T_cl + info['res'] / gain
            Tlist[bl][fq] = Tlist[bl].get(fq,[]) + [T * gain]
            Wlist[bl][fq] = Wlist[bl].get(fq,[]) + [gain]
            klist[bl] = ks
            #print fq, ks[fq][0]

    # Gotta do this one last time to catch the last integration.
    for sep,bls in sep2bl.items():
        if len(bls) < 2: continue
        for cnt,bl0 in enumerate(bls):
            if not Tlist.has_key(bl0) or len(Tlist[bl0]) == 0: continue
            #i0,j0 = a.miriad.bl2ij(bl0)
            #if 16 in [i0,j0]: print a.miriad.bl2ij(bl0), 'check'
            #print sep, [bl for bl in bls if len(Tlist[bl]) > 0]
            for bl1 in bls[cnt:]: # this includes "auto-pspecs"
                if len(Tlist[bl1]) == 0: continue
                for fq in Tlist[bl0]:
                    _Tlist = [Tlist[bl0][fq], Tlist[bl1][fq]]
                    _Wlist = [Wlist[bl0][fq], Wlist[bl1][fq]]
                    D, wgt = C.pspec.k3pk_from_Trms(_Tlist, _Wlist, k=klist[bl0][fq][0], B=B, fq=fq)
                    D,wgt = D[0], wgt[0]
                    D /= wgt # Write everything to file in mk^2 units
                    uvo['k3pk_fq'] = fq
                    uvo['k3pk_wgt'] = wgt
                    #uvo.copyvr(uvi)
                    print a.miriad.bl2ij(bl0), a.miriad.bl2ij(bl1), t, wgt
                    uvo.write((crd,curtime,(bl_index(bl0),bl_index(bl1))), D, n.zeros(D.shape, dtype=n.int))
    

