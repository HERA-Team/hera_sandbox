#! /usr/bin/env python

import aipy
import numpy
import capo
import os,sys
import optparse

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply_bandpass.py *uvcRREO')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

#Load passband
bp = numpy.load('bandpass.npz')['bandpass']

def mfunc(uv,p,d):
    d /= bp
    return p,d

#Apply passband
for file in args:
    print file, '->', file+'u'
    if os.path.exists(file+'u'):
        print '    File exists... skipping.'
        continue
    uvi = aipy.miriad.UV(file)
    uvo = aipy.miriad.UV(file+'u', status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,append2hist='OMNI_APPLY_PASSBAND:'+' '.join(sys.argv)+'\n')

