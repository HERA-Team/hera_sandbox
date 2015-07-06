#! /usr/bin/env python

import numpy
import aipy
import capo
import pylab
import sys

pol = 'xx'

aa = aipy.cal.get_aa('psa6622_v002',numpy.array([.15]))
sep, conj = capo.red.group_redundant_bls(aa.ant_layout)

bl_str = sep['0,2'] #change in row, change in col (ex: '0,1' is all the E/W 15 m bls)
goodant1 = 1
goodant2 = 4 #change these to reflect bl sep
bl_str = ['%d_%d' % aipy.miriad.bl2ij(bl) for bl in bl_str] #list comprehension
bls = ""
for i in bl_str:
    bls += i+','
bls = bls[:-1]
time, data, flags = capo.arp.get_dict_of_uv_data(sys.argv[1:],bls,pol)

bl_good = aipy.miriad.ij2bl(goodant1,goodant2)
ch = 100
good_data = data[bl_good][pol][:,ch]

for bl in data:
    comp_data = data[bl][pol][:,ch]
    if conj[bl]:
        comp_data = numpy.conj(comp_data)
 
    norm = numpy.sqrt(numpy.mean(numpy.abs(comp_data)**2)*numpy.mean(numpy.abs(good_data)**2))
    corr = numpy.mean(comp_data*numpy.conj(good_data))/norm #correlation
    corr_final = numpy.abs(corr)**2
    bl_final = aipy.miriad.bl2ij(bl)
    if corr_final < 0.9:
        print bl_final, corr_final

