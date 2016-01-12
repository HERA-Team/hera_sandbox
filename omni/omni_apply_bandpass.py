#! /usr/bin/env python

import aipy
import numpy
import capo
import os,sys
import optparse

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply_bandpass.py *uvcRRE.ucalbandpass.npz')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

def mfunc(uv,p,d):
    d /= bp
    return p,d

#Apply passband
for file in args:
    bp = numpy.load(file)['bandpass'] 
    files = numpy.load(file)['files']
    for fl in files:
        ifile = fl.split('.npz')[0]+'M' #applies bandpass to model vis UV file
        ofile = ifile+'u'
        print ifile, '->', ofile
        if os.path.exists(ofile):
            print '    File exists... skipping.'
            continue
        uvi = aipy.miriad.UV(ifile)
        uvo = aipy.miriad.UV(ofile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi,mfunc=mfunc,append2hist='OMNI_APPLY_PASSBAND:'+' '.join(sys.argv)+'\n')

