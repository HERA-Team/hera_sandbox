#! /usr/bin/env python

import aipy as a, sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
opts,args = o.parse_args(sys.argv[1:])

def mfunc(uv, p, d):
    return p, d
 
for filename in args:
    tfile = '.'.join([filename.split('.')[0]]+\
        [str(float('.'.join(filename.split('.')[1:3]))-2./24)] + \
        filename.split('.')[3:])+'C'#fixes timezone error in paper3
    print filename, '->', tfile
    uvi = a.miriad.UV(filename)
    if os.path.exists(tfile):
        print '    File exists... skipping.'
        continue
    uvo = a.miriad.UV(tfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,append2hist='Corrected time in filename')
