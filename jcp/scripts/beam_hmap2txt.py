#! /usr/bin/env python
import aipy as a, numpy as n

model = n.recfromtxt('/data1/paper/2010_beam/sdipole_05e_eg_ffx_150.txt',skiprows=2,usecols=[0,1,2])
th = model[:,0]*a.const.deg
ph = (model[:,1]-90)*a.const.deg

print min(th),max(th)
print min(ph),max(ph)

beam = a.map.Map(fromfits='rebeam.smoothe.fits')
val = 10*n.log(beam[th,ph])
val = n.where(val == -n.inf,0,val)

file = open('beammodel.txt','w')
for t,p,v in zip(th,ph,val):
    t = t / a.const.deg
    p = p / a.const.deg
    strng = str(t)+'           '+str(p)+'           '+str(v)+"\n"
    file.write(strng)

file.close()
    
