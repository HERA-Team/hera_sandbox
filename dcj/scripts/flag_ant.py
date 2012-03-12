#!/usr/bin/env python
#
#  apply_flag.py
#  
#
#  Created by Danny Jacobs on 6/16/10.
#  PAPER Project
#
"""
Flags any baselines, channels, polarizations selected at input.
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,os

o = optparse.OptionParser()
a.scripting.add_standard_options(o,ant=True,pol=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

for filename in args:
    print filename, '->', filename+'F'
    if os.path.exists(filename+'F'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
#    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
    flagants = map(a.miriad.bl2ij,[t[0] for t in a.scripting.parse_ants(opts.ant,uvi['nants'])])
#    chan = a.scripting.parse_chans(opts.chan,uvi['nchan'])
    uvo = a.miriad.UV(filename+'F', status='new')
    def mfunc(uv, p, d,f):
        if p[2] in flagants:
            return p,d,n.ones_like(f)
        else:
            return p,d,f
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc,raw=True,
        append2hist='FLAG ANT: Flagged baselines %s'%opts.ant)
    del(uvo)

        
