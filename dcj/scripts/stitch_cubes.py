#!/usr/bin/env python
#
#  stitch_cubes.py
#  
#
#  Created by Danny Jacobs on 9/28/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
from cosmo_units import *
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

root = '/data2/paper/sims/mcquinn/'
I6file = root+'zeta_evolve1_1000Mpc_400_lc_z6.1.ion.bool.npy'
I9file = root+'zeta_evolve1_1000Mpc_400_lc_z9.2.ion.bool.npy'
I14file = root+'zeta_evolve1_1000Mpc_400_lc_z14.7.ion.bool.npy'

I6 = n.load(I6file)

print I6.shape
I9 = n.load(I9file)
I14 = n.load(I14file)

#Ifile = 'zeta_evolve1_1000Mpc_800x2400_lc_z6_14.ion.bool.npy'
#I = n.load(Ifile)
#print I.shape
#X,Y,Z = n.mgrid[:400,:400,:200]
#X,Y,Z = n.mgrid[:800,:800,:400]
Z = n.zeros((400,400,200))
zi = n.arange(0,200)
for i in range(400):
    for j in range(400):
        Z[i,j,:] = zi
width = 200
wgt =  n.exp(-Z**2/(2*width**2))
#I = n.concatenate((I[:,:,:400], (I[:,:,400:800]*wgt + I[:,:,800:1200]*wgt[:,:,::-1])/2,
#        (I[:,:,1200:1600]*wgt + I[:,:,1600:2000]*wgt[:,:,::-1])/2,I[:,:,2000:]))
I = n.concatenate((I6[:,:,:200], (I6[:,:,200:]*wgt + I9[:,:,:200]*wgt[:,:,::-1])/2,
    (I9[:,:,200:]*wgt + I14[:,:,:200]*wgt[:,:,::-1])/2,I14[:,:,200:]),axis=2)
p.matshow(I[10,:,:].astype(n.bool))
p.show()