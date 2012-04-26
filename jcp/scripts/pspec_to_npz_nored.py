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

fq = .16
cen_fqs = n.array([fq])
B = .07
WINDOW = 'blackman-harris'
kwargs = {'cen_fqs':cen_fqs,'B':B, 'ntaps':1, 'window':WINDOW, 'bm_fqs':freqs.clip(.120,.190)}


#XXX these dictionaries will one day need polarization keys
Dat,Wgt = {},{}
tdat = {}
curtime, zen = None, None
for filename in args:
    print 'Reading', filename
    ofile = filename+'.nored.npz'
    if os.path.exists(ofile):
        print ofile, 'exists.  Skipping...'
        continue
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        bl = str(a.miriad.ij2bl(i,j))
        wbl = 'w'+str(a.miriad.ij2bl(i,j))
        #delay transform one sample at a time
        _d,ks =  C.pspec.Trms_vs_fq(freqs,d,**kwargs)
        _d,ks = _d[fq],ks[fq]
        try:
            #multiply adjacent time samples in delay space (leave out bias)
            Dat[bl] = Dat.get(bl,0) + tdat[bl]*_d.conj()
            Dat[wbl] = Dat.get(wbl,0) + 1.
        except(KeyError): pass
        tdat[bl] = _d

Dat['kpl'] = ks[1]

print 'Writing', ofile
n.savez(ofile,**Dat)
