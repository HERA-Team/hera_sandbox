#! /usr/bin/env python
"""
Select a range of jds from a list of input files
"""


import aipy as a, numpy as n, sys, os, optparse,re

o = optparse.OptionParser()
o.set_usage('jd_select.py [options] *.uv')
o.add_option('--jd_range',type=str,
     help='Range of dates in jd [ex: 2456789.1_2456999.2]')
o.add_option('--jdlist',type=str,
     help='a text file list of good jds.')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

#list of good jds.

if not opts.jd_range is None:
    jdmin,jdmax = map(float,opts.jd_range.split('_'))
    if jdmin>jdmax:
        raise(InputError, "invalid jd_range. {jd_range} jdmin>jdmax".format(jd_range=opts.jd_range))
else:
    raise(InputError,"please enter a jd_range")


def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])
JDs = [file2jd(os.path.basename(filename)) for filename in args]
for i in xrange(len(JDs)):
    if JDs[i]>jdmin and JDs[i]<jdmax: print args[i]

