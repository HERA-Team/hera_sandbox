#! /usr/bin/env python

import aipy, sys, os, numpy

skip_ants = [7, 13, 15]

def mfunc(uv, p, d, f):
    uvw, t, (i,j) = p
    if i in skip_ants or j in skip_ants: return p, None, None
    return p, d, f

for filename in sys.argv[1:]:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(filename+'c', status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='Removed antennas: %s\n' % (str(skip_ants)))
    del(uvo)

