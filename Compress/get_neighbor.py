#! /usr/bin/env python
"""
Returns a triplet of three adjacent files of sys.argv[1] in a space-separated string.
--- if there are no adjacent files, exit(1)

DFM
"""

from PDB import *
import sys
import re

fname=sys.argv[1]

pdb.verbose=False #just in case...
jd=re.findall(r'\d+\.\d+',fname)[0]

pol=fname.split('.')[-2]
try:
    assert(pol in ['xx','xy','yx','yy'])

    jd_lo=pdb.get('jd_lo','observations','JD',jd)[0][0]
    jd_hi=pdb.get('jd_hi','observations','JD',jd)[0][0]

    if jd_lo is None:
        print ""
    elif jd_hi is None:
        print ""
    else:
        lo_file=fname.replace(jd,jd_lo)
        hi_file=fname.replace(jd,jd_hi)
        print ' '.join([lo_file,fname,hi_file])
except(AssertionError):
    #this should happen for .npz files.
    print ""
