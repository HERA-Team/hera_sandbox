#!/usr/bin/env python

import aipy as a, sys, optparse, os
import numpy

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

def mfunc(uv,p,d,f):
    i,j = p[2]
    bl1 = str(i) + '_' + str(j)
    bl2 = str(j) + '_' + str(i)
    if bl1 in bls or bl2 in bls:
        return p,d,f
    else:
        return p,None,None

bls = []
for bl in opts.ant.split(','):
    bls.append(bl)

for filename in args:
    print filename, '->', filename+'A'
    uvi = a.miriad.UV(filename)
    if os.path.exists(filename+'A'):
        print '    File exists... skipping.'
        continue
    uvo = a.miriad.UV(filename+'A',status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,raw=True,mfunc=mfunc,append2hist='PULL ANTPOLS:'+' '.join(sys.argv)+'\n')


