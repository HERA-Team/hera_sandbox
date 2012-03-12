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
a.scripting.add_standard_options(o, ant=1,pol=1,chan=1)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

def gen_chans(chanopt, uv, coords, is_delay):
    """Return an array of active channels and whether or not a range of
    channels is selected (as opposed to one or more individual channels)
    based on command-line arguments."""
    is_chan_range = True
    if chanopt == 'all': chans = n.arange(uv['nchan'])
    else:
        chanopt = convert_arg_range(chanopt)
        if coords != 'index':
            if is_delay:
                def conv(c):
                    return int(n.round(c * uv['sdf'] * uv['nchan'])) \
                        + uv['nchan']/2
            else:
                def conv(c): return int(n.round((c - uv['sfreq']) / uv['sdf']))
        else:
            if is_delay:
                def conv(c): return int(c) + uv['nchan']/2
            else:
                def conv(c): return c
        chanopt = [map(conv, c) for c in chanopt]
        if len(chanopt[0]) != 1: 
            chanopt = [n.arange(x,y, dtype=n.int) for x,y in chanopt]
        else: is_chan_range = False
        chans = n.concatenate(chanopt)
    return chans.astype(n.int), is_chan_range

c1 = opts.chan.split(',')
chans = []
for c_spec in c1:
    cs = map(int,c_spec.split('_'))
    if len(cs)>1: chans.append(n.arange(cs[0],cs[-1]))
    else: chans.append(cs)
chans = n.concatenate(chans).astype(n.int32)
flags = n.ones_like(chans)

for filename in args:
    print filename, '->', filename+'r'
    if os.path.exists(filename+'r'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
#    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
#    flagants = map(a.miriad.bl2ij,[t[0] for t in a.scripting.parse_ants(opts.ant,uvi['nants'])])
#    chan = a.scripting.parse_chans(opts.chan,uvi['nchan'])
    uvo = a.miriad.UV(filename+'r', status='new')
    def mfunc(uv, p, d,f):
        f[chans] = n.logical_or(f[chans],flags)
        d[chans] = n.logical_not(flags)
        return p,d,f
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc,raw=True,
        append2hist='FLAG ANT: Flagged baselines %s'%opts.ant)
    del(uvo)