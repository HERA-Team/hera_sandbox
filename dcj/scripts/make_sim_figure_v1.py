#!/usr/bin/env python
#
#  make_sim_figure_v1.py
#  
#
#  Created by Danny Jacobs on 10/1/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
from cosmo_units import *

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

#make an inverse cosmology model
print "loading cosmology with O_l = %f, O_M = %f, Ho = %f"%(O_l,O_M,Ho)
zi = n.linspace(6,15)
D_i = [DM(z) for z in zi]
z_p = n.polyfit(D_i,zi,4)
z = n.poly1d(z_p)

DX = 1000 #distance in comoving mpc
print "loading ion"
I9 = n.fromfile('zeta_evolve1_1000Mpc_800_lc_z9.2.ion.bool',dtype=n.bool,count=800**2)
print "loading dens"
M9 = n.fromfile('zeta_evolve1_1000Mpc_800_lc_z9.2.dens',dtype='<f8',count=800**2)
I9.shape = M9.shape = (800,800)
C = lambda z: 26*n.sqrt((1+z)/10)

X,R = n.mgrid[:800,:800]
R = DM(9.2) + R*1000/800.

Z = z(R)
T = (1+M9)*I9*C(Z)
X *= 1000/800.

Theta = n.zeros((800,800))
for i in range(800):
    for j in range(800):
        Theta[i,j] = r2theta(X[i,j],Z[i,j]) * 180/n.pi #compute theta in degrees
p.axes(axisbg='k')
p.pcolor(Theta,Z,T,cmap='gist_heat',vmin=0)
p.colorbar(orientation='horizontal')
p.show()


