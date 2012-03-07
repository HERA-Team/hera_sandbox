#! /usr/bin/env python
import aipy as a, numpy as n
import sys, os, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True)
o.add_option('--gom', dest='gom', type='int',
    help='Index of gainometer antenna.')
o.add_option('-t', '--temp', dest='temp',
    help='Temperature data file')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)
temp=n.load(opts.temp)
T=temp['T']

for filename in args:
    dat = {}
    outfile = filename + 'T'
    print filename, '->', outfile
    if os.path.exists(outfile):
        print '    File exists: skipping'
        continue
    uv = a.miriad.UV(filename)
    uv.select('antennae', opts.gom, opts.gom)
    for (crd,t,(i,j)), d, f in uv.all(raw=True):
        factor = n.average(d.take(chans))
        dat[t] = factor
    del(uv)
    def mfunc(uv, p, d, f):
        crd,t,(i,j) = p
        if i == opts.gom or j == opts.gom: return p, None, None
        return p, d / dat[t], f

    uv = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uv)
    uvo.pipe(uv, mfunc=mfunc, raw=True,
        append2hist='GOM: gom=%d level=%f' % (opts.gom, opts.level))
    del(uvo)
