#! /usr/bin/env python
'''
Adds baselines together along columns, according to the redundant layout specified internally (PSA898).
'''

import aipy as a, numpy as n
import capo as C
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

# XXX Currently hardcoded for PSA898
A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]
ANTPOS = n.array([A_, B_, C_, D_])

bl2bin = {}
for i in range(ANTPOS.max()+1):
    for j in range(i, ANTPOS.max()+1):
        bl = a.miriad.ij2bl(i,j)
        offset = min(i % 4, j % 4) # for PSA898, this bins by column
        bin = a.miriad.ij2bl(i-offset, j-offset)
        bl2bin[bl] = bin

for filename in args:
    outfile = filename + 'S'
    print filename, '->', outfile
    if os.path.exists(outfile):
        print '    File exists, skipping.'
        continue
    print '    Summing baselines...'
    dsum,dwgt = {}, {}
    curtime = None
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        if t != curtime:
            curtime = t
            dsum[t], dwgt[t] = {}, {}
        pol = uvi['pol']
        if not dsum[t].has_key(pol):
            dsum[t][pol] = {}
            dwgt[t][pol] = {}
        bin = bl2bin[bl]
        #print (i,j),'->',a.miriad.bl2ij(bin)
        dsum[t][pol][bin] = dsum[t][pol].get(bin, 0) + n.where(f, 0, d)
        dwgt[t][pol][bin] = dwgt[t][pol].get(bin, 0) + n.logical_not(f).astype(n.int)
    uvi.rewind()

    print '    Writing output file'
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    def mfunc(uv, p, d, f):
        uvw,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        pol = uv['pol']
        try: _dsum,_dwgt = dsum[t][pol].pop(bl), dwgt[t][pol].pop(bl)
        except(KeyError): return p, None, None
        wgt = _dwgt.clip(1,n.Inf)
        f = n.where(_dwgt == 0, 1, 0)
        d = _dsum / wgt
        return (uvw,t,(i,j)), d, f

    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='PSPEC_COMBINE_BLS:' + ' '.join(sys.argv) + '\n')
    del(uvi); del(uvo)




