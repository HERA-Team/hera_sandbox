#!/usr/bin/env python
import aipy as a
import numpy as n
import sys,os

files = sys.argv[1:]

for f in files:
    startname = '.'.join(f.split('.')[1:3])
    if os.path.exists('%s_jds.txt'%startname):continue
    o = open('%s_jds.txt'%startname,'wb')
    uv = a.miriad.UV(f)
    uv.select('antennae',0,0)
    uv.select('polarization',-5,-5)
    for (crd,t,ij),d,f in uv.all(raw=True):
        o.write('%10.5f'%t+'\n')
    print 'writing %s'%startname
    o.close()
