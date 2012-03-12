#!/usr/bin/env python
#
#  blank_chans.py
#  
#
#  Created by Danny Jacobs on 5/25/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

for filename in args:
    print filename, '->', filename+'C'
    if os.path.exists(filename+'C'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    chan = a.scripting.parse_chans(opts.chan,uvi['nchan'])
    uvo = a.miriad.UV(filename+'C', status='new')
    curtime = 0
    def mfunc(uv, p, d,f):
        f[chan] |= 1
        return p,d,f    
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc,raw=True,
        append2hist='BLANK_CHANS: Flagged channels %s'%opts.chan)
    del(uvo)
