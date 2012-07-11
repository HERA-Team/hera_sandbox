#! /usr/bin/env python
import aipy as a, numpy as n
import capo as C
import sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, dec=True)
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    outfile = filename + 'C'
    print filename,'->',outfile
    if os.path.exists(outfile):
        print '    File exists.  Skipping...'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo._wrhd('history', uvi['history'] + 'PSPEC_SUM:' + ' ' .join(sys.argv) + '\n')
    
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
    uvi.select('decimate', opts.decimate, opts.decphs)

    dsum, dwgt = {},{}
    curtime = None
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        #if not '%d_%d' % (i,j) in opts.ant: continue # uv.select balks at large numbers of selections...
        if t != curtime:
            print curtime, dsum.keys()
            for fq in dsum:
                uvo['k3pk_fq'] = fq
                uvo['k3pk_wgt'] = dwgt[fq]
                d = dsum[fq] / dwgt[fq]
                uvo.write((crd,curtime,(0,0)), d, f)
            dsum, dwgt = {}, {}
            curtime = t
        fq = uvi['k3pk_fq']
        wgt = uvi['k3pk_wgt']
        dsum[fq] = dsum.get(fq,0) + d * wgt
        dwgt[fq] = dwgt.get(fq,0) + wgt

    # Gotta do this one more time to catch the last integration
    print curtime, dsum.keys()
    for fq in dsum:
        uvo['k3pk_fq'] = fq
        uvo['k3pk_wgt'] = dwgt[fq]
        d = dsum[fq] / dwgt[fq]
        uvo.write((crd,curtime,(0,0)), d, f)
