#! /usr/bin/env python
''' Sometimes files are not written to correctly. 
    Get segmentation faults, runtime errors, 
    Fatal Errors (Non-integral no. elemts in variable 
    npol, when scanning.
    This script checks uvfiles so that we dont get 
    these errors. 
'''

import aipy as a
import os, sys, optparse

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])
f = open('/data2/home/zakiali/psa_live/badpullantsfiles.txt','a')
for filename in args:
    print 'checking ', filename
    try: 
        uv = a.miriad.UV(filename)
    except Exception,e:
        print e, 'rerunning file'
        f.write('%s'%filename)
        os.system("rm -rf %s"%filename)
        

print 'done'




