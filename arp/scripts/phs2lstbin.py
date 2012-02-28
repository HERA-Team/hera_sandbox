#! /usr/bin/env python
"""
Rotate zenith UV data to a particular source.  Can specify 'zen' to phase data
to zenith, or nothing at all to just remove delay/offset phase components.
"""

import aipy as a, numpy as n, sys, os, optparse
import capo as C

o = optparse.OptionParser()
o.set_usage('phs2src.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

# A pipe to use for phasing to a source
curtime, zen = None, None
def phs(uv, p, d, f):
    global curtime
    uvw, t, (i,j) = p
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        u,v,lstbin = C.pspec.bin2uv(C.pspec.uv2bin(0,0,aa.sidereal_time()))
        zen = a.phs.RadioFixedBody(lstbin, aa.lat)
        zen.compute(aa)
    if i == j: return p, d, f
    d = aa.phs2src(d, zen, i, j)
    return p, d, f

# Process data
for filename in args:
    uvofile = filename + 'L'
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=phs, raw=True)

