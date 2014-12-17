#! /usr/bin/env python 
import numpy as n
import pylab as p
import sys, glob, os, optparse
import capo, aipy as a


args = sys.argv[1:]
aa = a.cal.get_aa('psa6240_v003', .1,.1,.1)
seps=['0,1','-1,1','1,1']
chan = [120] #channel to plot

sep2ij = {}
ij2sep = {}
s = ''
for sep in seps:
    sep2ij[sep] = capo.dfm.grid2ij(aa.ant_layout)[0][sep].split(',')
    s += capo.dfm.grid2ij(aa.ant_layout)[0][sep] + ','
    
    for ij in sep2ij[sep]:
        ij2sep[ij] = sep

conj = {}
toconj = capo.dfm.grid2ij(aa.ant_layout)[1]
for k in toconj.keys():
    conj['%d_%d'%a.miriad.bl2ij(k)] = toconj[k]



times,data,flags = capo.arp.get_dict_of_uv_data(files
