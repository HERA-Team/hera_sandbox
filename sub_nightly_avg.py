#! /usr/bin/env python
"""
Subtracts a nightly average to get to zero mean visibility
uv_avg > nightly_avg.py > sub_nightly_avg.py
"""
def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])

import aipy as a, numpy as n, sys, os, optparse, pickle,re
from smooth import smooth

o = optparse.OptionParser()
o.set_usage('uv_avg.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

for file in args:
    outfile = file +'x'
    print file, ' > ', outfile,
    if os.path.exists(outfile):
        print "File exists. Skipping."
    jd = file2jd(file)
    night = int(jd)
    night_avg_file = str(night)+'.avg.pkl'
    print night_avg_file
    if not os.path.exists(night_avg_file):
        print "please create nightly average: %s"%night_avg_file
        print "Skipping"
    uvi = a.miriad.UV(file)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    F = open(night_avg_file)
    AVG = pickle.load(F)
    def mfunc(uv,preamble,data):
        uvw,t,(i,j) = preamble
        bl = "%d_%d"%(i,j)
        data -= AVG.get(bl,0)
        return preamble, data
    uvo.pipe(uvi,mfunc=mfunc,
        append2hist='\n sub_nightly_avg.py: %s'%night_avg_file)
