#! /usr/bin/env python
"""
Select a range of integer jds (days) from a list of input files.
does a floor(jd) on the floating point filename jd. (bins files by day)
"""


import aipy as a, numpy as n, sys, os, optparse,re

o = optparse.OptionParser()
o.set_usage('jd_select.py [options] *.uv')
o.add_option('--jd_range',type=str,
     help='Range of dates in jd [ex: 2456789_2456999] ')
o.add_option('--jd_list',type=str,
     help='a text file list of good jds.')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])
if not opts.jd_list is None:
    jd_list = [int(l) for l in open(opts.jd_list).readlines() if not l.startswith('#')]

#list of good jds.    
if not opts.jd_range is None:
    jdmin,jdmax = map(float,opts.jd_range.split('_'))
    if jdmin>jdmax:
        raise(InputError, "invalid jd_range. {jd_range} jdmin>jdmax".format(jd_range=opts.jd_range))
    if not opts.jd_list is None:
        jd_list = [j for j in jd_list if (j<jdmax and j>jdmin)] #exclude anything outside our range
    else:
        jd_list = n.arange(jdmin,jdmax+1)
elif (not opts.jd_list is None) and (not opts.jd_range is None):
    raise(InputError,"please enter a jd_range and/or jd_list")


def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])
JDs = [file2jd(os.path.basename(filename)) for filename in args]
for i in xrange(len(JDs)):
    if n.floor(JDs[i]) in jd_list: print args[i]

