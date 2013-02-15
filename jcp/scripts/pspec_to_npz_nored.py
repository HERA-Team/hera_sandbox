#! /usr/bin/env python
import aipy as a, numpy as n, capo as C
import sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
freqs = uv['sfreq'] + n.arange(uv['nchan']) * uv['sdf']
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

fq = .1525
z = C.pspec.f2z(fq)
cen_fqs = n.array([fq])
B = .025
NEB = B / 2.006 #1-tap = 2.006; 3-tap = 1.756
#NEB = B / 1.756 #1-tap = 2.006; 3-tap = 1.756
WINDOW = 'blackman-harris'
kwargs = {'cen_fqs':cen_fqs,'B':B, 'ntaps':1, 'window':WINDOW, 'bm_fqs':freqs.clip(.120,.190)} 
bm = n.polyval(C.pspec.DEFAULT_BEAM_POLY,fq)

#XXX these dictionaries will one day need polarization keys
Dat,Wgt = {},{}
tdat = {}
curtime, zen = None, None
for filename in args:
    print 'Reading', filename
    ofile = filename+'.%s.nored.npz' % str(B)
    if os.path.exists(ofile):
        print ofile, 'exists.  Skipping...'
        continue
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            old_zen = zen
            aa.set_jultime(t)
            lst = aa.sidereal_time()
            zen = a.phs.RadioFixedBody(lst, aa.lat)
            zen.compute(aa)
            curtime = t
        bl = str(a.miriad.ij2bl(i,j))
        wbl = 'w'+str(a.miriad.ij2bl(i,j))
        #flip samples across the uv plane
        u,v,w = aa.gen_uvw(i,j,src=zen)
        u = u.flatten()[-1]
        conj = False
        if u < 0: conj = True
        #phase second integration to first
        td = d.copy()
        #if conj: n.conj(td)
        _td,tks = C.pspec.Trms_vs_fq(freqs,td,**kwargs)
        _td,tks = _td[fq],tks[fq]
        if old_zen != None:
            d *= n.conj(aa.gen_phs(zen,i,j)) * aa.gen_phs(old_zen,i,j)
        #if conj: d = n.conj(d)
        #delay transform one sample at a time
        _d,ks =  C.pspec.Trms_vs_fq(freqs,d,**kwargs)
        _d,ks = _d[fq],ks[fq]
        scalar = C.pspec.X2Y(z) * bm * NEB
        try:
            #multiply adjacent time samples in delay space (leave out bias)
            Dat[bl] = Dat.get(bl,0) + scalar * tdat[bl]*_d.conj()
            Dat[wbl] = Dat.get(wbl,0) + 1.
        except(KeyError): pass
        tdat[bl] = _td

    Dat['kpl'] = ks[1]

    print 'Writing', ofile
    n.savez(ofile,**Dat)
