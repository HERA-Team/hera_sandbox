#! /usr/bin/env python

import aipy as a, numpy as n, sys, os, ephem, optparse, glob
from time import strftime, gmtime, mktime, strptime

o = optparse.OptionParser()
o.set_usage('miriad_prep.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
for filename in args:
    print filename,'->',filename+'M'
    if os.path.exists(filename+'M'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'M', status='new')
    antpos = n.array([ant.pos for ant in aa])
    antpos.shape  = antpos.shape[0]*antpos.shape[1]
    uvo.add_var('restfreq','d')
    uvo['restfreq'] = uv['sfreq']
    override = {'antpos':antpos}
    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi,append2hist='miriad_prep.py: Added correct antenna positions for MIRIAD export.')
