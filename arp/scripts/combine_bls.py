#! /usr/bin/env python
import aipy as a, numpy as n
import capo as C
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
ants = a.scripting.parse_ants(opts.ant, uv['nants'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

UV_RES = 2.
LST_RES = 2*n.pi*a.ephem.second

for filename in args:
    outfile = filename + 'S'
    print filename, '->', outfile
    if os.path.exists(outfile):
        print '    File exists, skipping.'
        continue
    dsum,dwgt = {}, {}
    bin2bl = {}
    curtime = None
    print '    Summing baselines...'
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, ants, opts.pol)
    #uvi.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        if t != curtime:
            curtime = t
            aa.set_jultime(t)
        u,v,w = aa.gen_uvw(i, j, 'z')
        u,v = u[0,-1],v[0,-1]
        if u < -UV_RES/2: # Conjugate to fold UV plane to u > 0
            u,v,d = -u,-v,n.conj(d)
        # XXX should bin then unbin to determine folding (and flip -v for u=0)
        bin = C.pspec.uv2bin(u,v, aa.sidereal_time(), uv_res=UV_RES, lst_res=LST_RES)
        bin2bl[bin] = bin2bl.get(bin,[]) + [bl]
        dsum[bin] = dsum.get(bin, 0) + n.where(f, 0, d)
        dwgt[bin] = dwgt.get(bin, 0) + n.logical_not(f).astype(n.int)
    uvi.rewind()
    for bin in bin2bl.keys():
        bin2bl[bin].sort()
        bin2bl[bin] = bin2bl[bin][0]

    print '    Writing output file'
    curtime = None
    def mfunc(uv, p, d, f):
        global curtime
        uvw,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        if t != curtime:
            curtime = t
            aa.set_jultime(t)
        u,v,w = aa.gen_uvw(i, j, 'z')
        u,v = u[0,-1],v[0,-1]
        if u < -UV_RES/2: u,v = -u,-v
        bin = C.pspec.uv2bin(u,v, aa.sidereal_time(), uv_res=UV_RES, lst_res=LST_RES)
        if bl != bin2bl[bin]: return (uvw,t,(1,1)), None, None
        wgt = dwgt[bin].clip(1,n.Inf)
        f = n.where(dwgt[bin] == 0, 1, 0)
        d = dsum[bin] / wgt
        return (uvw,t,(i,j)), d, f

    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='COMBINE BLS:' + ' '.join(sys.argv) + '\n')
    del(uvi); del(uvo)




