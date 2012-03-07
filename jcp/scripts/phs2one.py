#! /usr/bin/env python

import aipy as a, sys, optparse, os, numpy as np

o = optparse.OptionParser()
a.scripting.add_standard_options(o, src=True, cal=True, ant=True)
o.add_option('--conj', dest='conj', action='store_true',
    help='Conjugate all the data before phasing')         
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)

nbls = a.scripting.parse_ants(opts.ant, uv['nants'])
if len(nbls) != 1: raise ValueError("select only one baseline!")
bl = nbls[0][0]
ants =  a.miriad.bl2ij(bl)
l,m = ants

def mfunc(uv, p, d):
    crd, t, (i, j) = p
    aa.set_jultime(t)
    cat.compute(aa)
    if opts.conj:
        d = np.conj(d)
    d = d * aa.gen_phs(cat[k], l, m)
    p = crd, t, (i, j)
    return p, d


for k in cat:
    for filename in args:
        outfile = filename+'.phs%d%d'%(l,m)+k
        if opts.conj: outfile = outfile + '.conj'
        print filename, '->', outfile
        uvi = a.miriad.UV(filename)
        if os.path.exists(outfile):
            print '    File exists... skipping.'
            continue
        uvo = a.miriad.UV(outfile,status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi,mfunc=mfunc,append2hist='Phased to antenna %d%d.\n' % (l,m))
