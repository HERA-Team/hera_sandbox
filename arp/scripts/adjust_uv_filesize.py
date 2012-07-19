#! /usr/bin/env python
"""
"""

import aipy as a, sys, os, numpy as n, optparse

o = optparse.OptionParser()
o.set_usage('adjust_uv_filesize.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-n', '--nint', type='int',
    help='Number of contiguous integrations to output per file.')
opts, args = o.parse_args(sys.argv[1:])

uvo = None
times = []
for uvfile in args:
    fileparts = uvfile.split('.')
    filename_fmt = '.'.join([fileparts[0], '%13.5f', fileparts[-1]+'z'])
    uvi = a.miriad.UV(uvfile)
    if not uvo is None: print uvfile,'->',outfile
    
    for p,d,f in uvi.all(raw=True):
        t = p[1]
        # Check that deltas between jds are consistent and end file otherwise
        if len(times) % (opts.nint-1) == 0 or (times[-1] != t and len(times) > 1 \
                and n.abs((t-times[-1])/(times[1]-times[0])-1) > 0.5):
            if len(times) > 1 and times[-1] != t and n.abs((t-times[-1])/(times[1]-times[0])-1) > 0.5:
                print '    Time gap detected.  Ending current file.'
            if not uvo is None: del(uvo) # Close out existing file
            outfile = filename_fmt % t
            print uvfile,'->',outfile
            uvo = a.miriad.UV(outfile, status='new')
            uvo.init_from_uv(uvi)
            uvo._wrhd('history', uvo['history'] + ' '.join(sys.argv) + '\n')
            times = []
        if len(times) == 0 or times[-1] != t: times.append(t)
        uvo.copyvr(uvi)
        uvo.write(p,d,f)
del(uvo) 
