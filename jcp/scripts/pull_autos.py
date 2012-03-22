#! /usr/bin/env python

import aipy as a, sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    print filename, '->', filename+'A'
    uvi = a.miriad.UV(filename)
    if os.path.exists(filename+'A'):
        print '    File exists... skipping.'
        continue
    a.scripting.uv_selector(uvi, ants=opts.ant)
    uvo = a.miriad.UV(filename+'A',status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,append2hist='Selected autocorrelations\n')
