#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse

model = n.recfromtxt('/data1/paper/2010_beam/sdipole_05e_eg_ffx_150.txt',skiprows=2,usecols=[0,1,2])
#model = n.recfromtxt('beammodel.txt',usecols=[0,1,2])
th = model[:,0]*a.const.deg
ph = (model[:,1])*a.const.deg
val = model[:,2]
max = n.max(10**(val/10))
print n.min(th),n.max(th)
print n.min(ph),n.max(ph)

beam = a.map.Map(nside=64,interp=True)
beam.add((th,ph),1,(10**(val/10)/max))
#beam.add((th,ph),1,val)
#beam.add((th,ph),1,10*n.log10(val))
beam.to_fits('model150mhz_n64.fits',clobber=True)
