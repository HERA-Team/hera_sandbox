"""
casapy script for checking measurement sets and re-importing them from uvfits if necessary
"""

from time import time
import numpy as n
import aipy as a
import sys,os,optparse
import re,shutil
t0 = time()

import mwapy
import mwapy.get_observation_info
from mwapy.obssched.base import schedule

#db=schedule.getdb()

####################################
##     Parse inputs               ##
####################################
o = optparse.OptionParser()
o.set_usage('fix_ms.py [options] <obsid>')
o.set_description(__doc__)
o.add_option('--restore_flags',type=str,default='MWAflag',
    help='Name of flag table to restore from backup should the ms be corrupted.')
for i in range(len(sys.argv)):
    if sys.argv[i]==inspect.getfile( inspect.currentframe()):break
opts, args = o.parse_args(sys.argv[i+1:])

for obsid in args:
    uvfits = obsid+'/'+obsid+'.uvfits'
    vis = obsid +'/'+obsid+'.ms'
    print "testing ",vis,
    try: 
        ms.open(vis)
        print "[PASS]"
    except(StandardError):
        print "[FAIL]",
        print "checking for existing ms",
        if os.path.exists(vis):
            print "file exists but must be bad, deleting",
            shutil.rmtree(vis)
        print "re-importing from ",uvfits,
        importuvfits(fitsfile=uvfits,vis=vis)
        if not opts.restore_flags is None:
            print "restoring flag version, ",opts.restore_flags
            flagmanager(vis=vis,mode='restore',versionname=opts.restore_flags)
        else: print 


