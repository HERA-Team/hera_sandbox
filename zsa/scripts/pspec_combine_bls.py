#! /usr/bin/env python
'''
Adds baselines together along columns, according to the redundant layout specified internally (PSA898).
Updating script to PSA6240, which is not in bit reversed order.
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
#E_ = [i+4 for i in A_]
#F_ = [i+5 for i in A_]
#G_ = [i+6 for i in A_]
#H_ = [i+7 for i in A_]
ANTPOS = n.array([A_, B_, C_, D_])

# PSA-128, JD2456240...
A_ = [49,41,47,19,29,28,34,51]
B_ = [10, 3,25,48,24,55,27,57]
C_ = [ 9,58, 1, 4,17,13,56,59]
D_ = [22,61,35,18, 5,32,30,23]
E_ = [20,63,42,37,40,14,54,50]
F_ = [43, 2,33, 6,52, 7,12,38]
G_ = [53,21,15,16,62,44, 0,26]
H_ = [31,45, 8,11,36,60,39,46]
ANTPOS_6240 = n.array([A_, B_, C_, D_,E_,F_,G_,H_])

ANTPOS = ANTPOS_6240

bl2bin = {}
conjbl = {}
for x in range(ANTPOS.size):
    for y in range(x, ANTPOS.size):
        i = ANTPOS.flat[x]
        j = ANTPOS.flat[y]
        bl = a.miriad.ij2bl(i,j)
        #this bins by column
        offset = min(x % ANTPOS.shape[0], y % ANTPOS.shape[0]) 
        ox, oy = x-offset, y-offset
        #match baseline orientation of bin
        oi, oj = ANTPOS.flat[ox], ANTPOS.flat[oy]
        bin = a.miriad.ij2bl(oi, oj)
        conjbl[bl] = ((i < j) != (oi < oj)) #check match bl orientation
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
        if conjbl[bl]: d = d.conj()
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
        #bl is equal to bin only for the top bl in a column.
        try: _dsum,_dwgt = dsum[t][pol].pop(bl), dwgt[t][pol].pop(bl)
        except(KeyError): return p, None, None
        wgt = _dwgt.clip(1,n.Inf)
        f = n.where(_dwgt == 0, 1, 0)
        d = _dsum / wgt
        return (uvw,t,(i,j)), d, f

    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='PSPEC_COMBINE_BLS:' + ' '.join(sys.argv) + '\n')
    del(uvi); del(uvo)




