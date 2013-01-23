#! /usr/bin/env python
import aipy as a, numpy as n
import capo as C
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
ants = a.scripting.parse_ants(opts.ant, uv['nants'])
#aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

for filename in args:
    outfile = filename + 'P'
    print filename, '->', outfile
    if os.path.exists(outfile):
        print '    File exists, skipping.'
        continue
    dsum,dwgt = {}, {}
    curtime = None
    print '    Summing %s...' % opts.pol
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, ants, opts.pol)
    #uvi.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        if t != curtime:
            curtime = t
            #aa.set_jultime(t)
            dsum[t], dwgt[t] = {}, {}
        dsum[t][bl] = dsum[t].get(bl, 0) + n.where(f, 0, d)
        dwgt[t][bl] = dwgt[t].get(bl, 0) + n.logical_not(f).astype(n.int)
    uvi.rewind()

    print '    Writing output file'
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    def mfunc(uv, p, d, f):
        uvw,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        try: _dsum,_dwgt = dsum[t].pop(bl), dwgt[t].pop(bl)
        except(KeyError): return p, None, None
        #uvo['pol'] = a.miriad.str2pol['I'] # this doesn't work, so don't bother
        wgt = _dwgt.clip(1,n.Inf)
        f = n.where(_dwgt == 0, 1, 0)
        d = _dsum / wgt
        return (uvw,t,(i,j)), d, f

    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='COMBINE_POL:' + ' '.join(sys.argv) + '\n')
    del(uvi); del(uvo)




