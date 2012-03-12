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
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=1,pol=1,chan=1)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

chans = 
def flag_pipe(uv,p,d,f):
    
    

for file in args:
    uvofile = uvfile+'F'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)
    chans = a.scripting.parse_chans(opts.chan, uvi['nchan'])
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for p,d,f in uv.all(raw=True):
        f |= n.ones_like(f)
        