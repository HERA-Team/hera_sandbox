#! /usr/bin/env python

import aipy
import numpy
import capo
import os,sys
import optparse

### Options ###
o = optparse.OptionParser()
o.set_usage('apply_uCal.py *uCalResults.npz *uvcRRE')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

def mfunc(uv,p,d,f):
    d /= bp
    return p,d,f

npz = numpy.load(args[0])
bp = npz['bandpassFit']

#Apply bandpass
for file in args[1:]:
    ofile = file+'u'
    print file, '->', ofile
    if os.path.exists(ofile):
        print '    File exists... skipping.'
        continue
    uvi = aipy.miriad.UV(file)
    uvo = aipy.miriad.UV(ofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,raw=True,mfunc=mfunc,append2hist='APPLY_UCAL:'+' '.join(sys.argv)+'\n')


