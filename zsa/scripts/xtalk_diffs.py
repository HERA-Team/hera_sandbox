#!/usr/bin/env python
import aipy as a
import numpy as n 
import pylab as p
import sys,optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
o.add_option('--sep', action='store', default='1,0',
              help='Which separation to look at')
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa('psa6240_v003', n.array([.150]))

ANTPOS = aa.ant_layout
#make dict of separations.
bl2sep = {}
sep2bl = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
                if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                #sep = a.miriad.ij2bl(rj-ri, cj-ci)
                sep = (cj-ci, rj-ri) #(dx,dy) in row spacing units
                i,j = ANTPOS[ri,ci], ANTPOS[rj,cj]
                bl = a.miriad.ij2bl(i,j)
                if i > j: 
                    i,j = j,i
                    sep = (sep[0]*-1,sep[1]*-1)
                bl2sep[bl] = sep
                sep2bl[sep] = sep2bl.get(sep,[]) + [bl]

mysep = tuple([int(i) for i in opts.sep.split(',')])
print 'looking at baselines with separation ', mysep
print 'baselines are : ', sep2bl[mysep], 'len = ', len(sep2bl[mysep])
#turn the list in to a string that aipy likes to choose ants with
s = ''.join([str(a.miriad.bl2ij(bl)) for bl in sep2bl[mysep]])
s = s.replace(',', '_')
s = s.replace('(','')
s = s.replace(')', ',')
s = s.replace(' ', '')
print s

#the diffs of different baselines, based on separation types.
diffs = []
all_bls = [] #holds data for all the bls for a sep, for a time slice.
tslice=0.
for file in args:
    print file
    uv = a.miriad.UV(file)
    a.scripting.uv_selector(uv, ants=s, pol_str='xx')
    for (uvw,t,(i,j)),d in uv.all():
        if (t != tslice) and (len(all_bls) > 0): 
            all_bls = n.ma.array(all_bls)  
            all_bls1 = all_bls[::2]
            all_bls2 = all_bls[1::2]
            diff = all_bls2 - all_bls1
            diffs.append(n.ma.average(diff, axis=0))
            tslice = t
            all_bls = []
        else:
            all_bls.append(d)

print diffs
p.plot(diffs)
