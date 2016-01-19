#! /usr/bin/env python
#plot the selected spacings in the input omnical solutions
#
# usage: 

from capo.omni import from_npz
from capo.dfm import grid2ij
import sys,os
from pylab import *
import numpy as n
import ipdb as PDB
cal = 'psa6622_v002'
exec('from {cal} import prms'.format(cal=cal))
grid=grid2ij(prms['ant_layout'])[0]
def ij2grid(bl):
    for sep in grid:
        if bl in grid[sep].split(','):
            return sep
def bltxt(bltuple):
    return '%s_%s'%bltuple
print ij2grid('0_103')
seps = ['0,1']
data = {}
for filename in sys.argv[1:]:
    meta, gains, vismdl, xtalk, jds, lsts, freqs = from_npz(filename)
    #for select off my bls
    for pol in vismdl:
        try: data[pol]
        except(KeyError): data[pol] = {} 
        for bl in vismdl[pol]:
            blsep = ij2grid(bltxt(bl))
            if blsep in seps:
                try: data[pol][blsep].append(vismdl[pol][bl])
                except(KeyError):data[pol][blsep] = [vismdl[pol][bl]]
datacount = 0
for pol in data:
    for sep in data[pol]:
        try:
            data[pol][sep] = n.concatenate(data[pol][sep])
        except(TypeError):
            PDB.set_trace()
        datacount += 1
figure()
rows = n.ceil(n.sqrt(datacount))
cols = n.ceil(datacount/rows)
i=1
print data
for pol in data:
    for sep in data[pol]:
        subplot(rows,cols,i)
        imshow(n.angle(data[pol][sep]),aspect='auto')
        title(sep+' '+pol)
        i+=1

show()
