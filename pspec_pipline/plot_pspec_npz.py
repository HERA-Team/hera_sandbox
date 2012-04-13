#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys

PLOT_SPEC = False
PLOT_PSPEC = False
B = .008
WINDOW = 'blackman-harris'

#if True:
#    sdf = freqs[1] - freqs[0]
#    add_chan = 1024
#    lo_fq = n.arange(freqs[0] - add_chan*sdf,freqs[0],sdf)
#    hi_fq = n.arange(freqs[-1]+sdf,freqs[-1] + add_chan*sdf,sdf)
#    freqs = n.concatenate([lo_fq, freqs, hi_fq])
#    sampling = n.concatenate([n.zeros_like(lo_fq), sampling, n.zeros_like(hi_fq)])
#    window = n.concatenate([n.zeros_like(lo_fq), window, n.zeros_like(hi_fq)])
#
#_window = n.fft.ifft(window * sampling)

fdat,fwgt = {}, {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    f = n.load(filename)
    freqs = f['freqs']
    bins = f['bins']
    bmax, bwgt = None, 0.
    for cnt,b in enumerate(f['bins']):
        wgt = n.sum(f['wgt'][cnt], axis=0).max()
        if wgt > bwgt: bmax,bwgt = b,wgt
        fdat[b] = fdat.get(b,0) + f['dat'][cnt]
        fwgt[b] = fwgt.get(b,0) + f['wgt'][cnt]
    print '   ', len(bins), bmax, bwgt
print len(fdat)

Dsum,Dwgt = {}, {}
kwargs = {'cen_fqs':n.arange(.115,.190,.005),'B':B, 'ntaps':3, 'window':WINDOW, 'bm_fqs':freqs.clip(.120,.190)}
window = a.dsp.gen_window(freqs.size, window=WINDOW)
for cnt,b in enumerate(fdat):
    ubin,vbin,lstbin = C.pspec.bin2uv(b)
    umag = n.sqrt(ubin**2 + vbin**2)
    umag = 2**int(n.around(n.log2(umag.clip(0.5,n.Inf))))
    #if umag != 10: continue
    d, w = n.sum(fdat[b], axis=0), n.sum(fwgt[b], axis=0)
    wgt = w.max()
    print cnt, b, (ubin,vbin), umag, wgt
    #if wgt < 500: continue
    #if wgt < 200: continue
    if wgt < 100: continue
    #if wgt < 50: continue

    if PLOT_SPEC:
        d /= wgt; w /= wgt
        #_d = n.fft.ifft(d*window)
        #_w = n.fft.ifft(w*window)
        #_d_cl, info = a.deconv.clean(_d, _w, tol=1e-9, stop_if_div=False)
        #d_mdl = n.fft.fft(_d_cl) 
        #Tmdl,ks = C.pspec.Trms_vs_fq(freqs, d_mdl, **kwargs)
        jy2T = C.pspec.jy2T(freqs.clip(.120,.190))
        p.plot(freqs, d.real*jy2T)
        #p.plot(freqs, d_mdl.real*jy2T)
        #p.plot(freqs, (d-d_mdl).real*jy2T)
        p.show()
    
    Tlist,Wlist = {}, {}
    for d,w in zip(fdat[b], fwgt[b]):
        try: wgt = w.max()
        except: continue
        d /= wgt; w /= wgt
        if False: d_res = d - d_mdl * w
        else: d_res = d
        Tres,ks = C.pspec.Trms_vs_fq(freqs, d_res, **kwargs)
        W,ks = C.pspec.Trms_vs_fq(freqs, w, **kwargs)
        
        for fq in Tres:
            gain = n.abs(W[fq][0])
            T_cl, info = a.deconv.clean(Tres[fq], W[fq], tol=1e-9, maxiter=100,
                stop_if_div=False, verbose=False)
            #print '   ', int(n.around(fq*1e3)), info['term']
            if False: T = Tmdl[fq] + T_cl + info['res']
            else: T = T_cl + info['res'] / gain
            Tlist[fq] = Tlist.get(fq,[]) + [T * wgt]
            Wlist[fq] = Wlist.get(fq,[]) + [wgt]
    for fq in Tlist:
        if not Dsum.has_key(fq): Dsum[fq], Dwgt[fq] = {}, {}
        if len(Tlist[fq]) < 2: continue
        D, wgt = C.pspec.k3pk_from_Trms(Tlist[fq], Wlist[fq], k=ks[fq][0], B=B, fq=fq)
        if PLOT_PSPEC and not PLOT_SPEC:
            p.loglog(ks[fq][0], n.abs(D.real).clip(1e0,n.Inf), ',', label=str(b))
        print '   ', fq, 'adding to D^2 with weight:', wgt
        Dsum[fq][umag] = Dsum[fq].get(umag,0) + D * wgt
        Dwgt[fq][umag] = Dwgt[fq].get(umag,0) + wgt

for fq in Dsum:
    k = ks[fq][1]
    D = {}
    D['k'] = k
    for umag in Dsum[fq]:
        print '   ', int(1e3*fq), umag
        D[str(umag)] = (Dsum[fq][umag] / Dwgt[fq][umag])
    n.savez('B_%s.npz' % fq, **D)
        #p.loglog(k, n.abs(D.real).clip(1e0,n.Inf), 
        #    label='%d,%d' % (int(1e3*fq), umag))
        #p.loglog(k, n.abs(D.imag).clip(1e0,n.Inf), 
        #    label='%d,%d' % (int(1e3*fq), umag))

p.legend(loc='upper left')
p.show()
        
