#! /usr/bin/env python
"""
Averages the outputs of uv_avg into a single nightly file
"""


import aipy as a, numpy as n, sys, os, optparse, pickle,re
from smooth import smooth

o = optparse.OptionParser()
o.set_usage('nightly_avg.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])


#break the list of files into nights
def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])
nights = {}
nights_c  = {}
#nights = set([round(file2jd(file),0) for file in args])
#nights = dict(zip(nights,[[]]*len(nights)))
for file in args:
    print file
    jd = file2jd(file)
    night = int(jd)
    F = open(file)
    AVG = pickle.load(F)
    if not nights.has_key(night): 
        nights[night] = {}
        nights_c[night] = {}
    for bl in AVG:
        if bl=='freqs':freqs = AVG['freqs'];continue
        nights[night][bl] = nights[night].get(bl,0) + AVG[bl]
        nights_c[night][bl] = nights_c[night].get(bl,0) + 1
for night in nights:
    for bl in nights[night]:
        N = nights_c[night][bl]
        if N>0: nights[night][bl] /= N
    nights[night]['freqs'] = freqs
for night in nights:
    outfile = str(night)+'.avg.pkl'
    print "writing: ",outfile
    F = open(outfile,'w')
    pickle.dump(nights[night],F)
    F.close()
