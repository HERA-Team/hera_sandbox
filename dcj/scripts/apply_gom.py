#!/usr/global/paper/bin/python
"""
Divide out the gom.
"""

__version__ = '0.0.1'

import aipy as a, numpy as n, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('apply_gom.py [options] *.uv')
o.set_description(__doc__)
opts, args = o.parse_args(sys.argv[1:])

def mdl(uv, p, d, f):
    
# Process all files passed from the command line.
for filename in args:
    print filename,'->',filename+'g'
    if os.path.exists(filename+'b'):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'b', status='new')
    uvo.init_from_uv(uvi, exclude=ignore_vars)
    for p,d in uvi.all():
        
    nchan = uvi['nchan']
    nants = uvi['nants']
    try:
        # If there is a bandpass item, we're going apply it to data
        bp = uvi['bandpass'].real   # Sync'd with pocket_corr.py in corr pkg.
        bp.shape = (nants, nchan)
        print
    except:
        print 'No bandpass found'
        bp = n.ones((nants, nchan))
        print '.'
    def f(uv, preamble, data, flags):
        uvw, t, (i,j) = preamble
        if i == j: data = n.polyval(cpoly, data)
        data *= bp[i,:] * bp[j,:] * opts.scale
        return preamble, data, flags
    uvo.pipe(uvi, mfunc=f, raw=True,
        append2hist='APPLY_BP: ver=%s, corr type=%s, scale=%f\n' % \
            (__version__, opts.linearization, opts.scale))
