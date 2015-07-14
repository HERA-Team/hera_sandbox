#! /usr/bin/env python

import aipy as a, numpy as n, sys, os, optparse,re

o = optparse.OptionParser()
o.set_usage('list_jds.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])
def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])

JDs = [file2jd(os.path.basename(filename)) for filename in args]
JD_ints = list(set(map(n.floor,JDs)))
for jd in JD_ints:
    print int(jd),
