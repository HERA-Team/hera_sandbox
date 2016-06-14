#! /usr/bin/env python
"""
Plots the output of  uv_avg in a waterfall
"""


import aipy as a, numpy as n, sys, os, optparse, pickle,re
from smooth import smooth
from pylab import *
o = optparse.OptionParser()
o.set_usage('plot_xtalk_waterfall.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])
def file2jd(file):
    return float(re.findall('\D+(\d+.\d+)\D+',file)[0])

mybl= '41_49'
waterfall = []
times = []
freqs = []
for filename in args:
    print filename
    AVG = pickle.load(open(filename))
    waterfall.append(AVG[mybl]/AVG['counts'][mybl].astype(float))
    times.append([file2jd(filename)]*len(AVG['freqs']))
    freqs.append(AVG['freqs'])
T,F,D = n.array(times),n.array(freqs),n.array(waterfall)
print T.shape,F.shape,D.shape
pcolor(T,F,D)
n.savez('xtalk_waterfall.npz',T=T,F=F,D=D)
show()


##break the list of files into nights
#nights = {}
#nights_c  = {}
##nights = set([round(file2jd(file),0) for file in args])
##nights = dict(zip(nights,[[]]*len(nights)))
#for file in args:
#    print file
#    jd = file2jd(file)
#    night = int(jd)
#    F = open(file)
#    AVG = pickle.load(F)
#    if not nights.has_key(night): 
#        nights[night] = {}
#        nights_c[night] = {}
#    for bl in AVG:
#        if bl=='freqs':freqs = AVG['freqs'];continue
#        if bl=='counts':continue
#        nights[night][bl] = nights[night].get(bl,0) + AVG[bl]*AVG['counts'][bl].astype(float)
#        nights_c[night][bl] = nights_c[night].get(bl,0) + AVG['counts'][bl]
#for night in nights:
#    for bl in nights[night]:
#        N = nights_c[night][bl]
#        nights[night][bl][N>0] /= N[N>0]
#    nights[night]['counts'] = nights_c[night]
#    nights[night]['freqs'] = freqs
#for night in nights:
#    outfile = str(night)+'.avg.pkl'
#    print "writing: ",outfile
#    F = open(outfile,'w')
#    pickle.dump(nights[night],F)
#    F.close()
