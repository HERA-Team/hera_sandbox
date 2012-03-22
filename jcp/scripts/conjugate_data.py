#! /usr/bin/env python

import aipy as a, sys, optparse, os, numpy as np

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
opts,args = o.parse_args(sys.argv[1:])

def mfunc(uv, p, d):
    d = np.conj(d)
    return p, d
 
for filename in args:
    print filename, '->', filename+'.conj'
    uvi = a.miriad.UV(filename)
    if os.path.exists(filename+'.conj'):
        print '    File exists... skipping.'
        continue
    uvo = a.miriad.UV(filename+'.conj',status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,append2hist='Conjugated all the data.\n')
