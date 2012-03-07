#! /usr/bin/env python

import aipy as a, sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
o.add_option('-t', '--time', dest='time', default=0, type='float',
    help='Time correction value in hours.  Negative to subtract.')
opts,args = o.parse_args(sys.argv[1:])

def mfunc(uv, p, d):
    crd, t, (i, j) = p
    t += opts.time/24
    p = crd, t, (i, j)
    return p, d
 
for filename in args:
    print filename, '->', filename+'T'
    uvi = a.miriad.UV(filename)
    if os.path.exists(filename+'T'):
        print '    File exists... skipping.'
        continue
    uvo = a.miriad.UV(filename+'T',status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,append2hist='Corrected timestamps by %d \n' % opts.time)
