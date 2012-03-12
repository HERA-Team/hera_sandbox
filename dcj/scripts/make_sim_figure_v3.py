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
import matplotlib as mpl
#mpl.rcParams['text.color'] = 'w'
#mpl.rcParams['ytick.color'] = 'w'
#mpl.rcParams['xtick.color'] = 'w'
#mpl.rcParams['axes.edgecolor'] = 'w'
#mpl.rcParams['savefig.facecolor'] = 'w'
#mpl.rcParams['savefig.edgecolor'] = 'w'

mpl.rcParams['font.weight'] = 'medium'
"""
I put two files in your directory:  1024_1200Mpc_lc_z7.1.ion and the same but 
for density (it's currently being transferred).  Could you look at this and tell 
me if it looks okay?  It should be one 1200 Mpc box that contains most of the
 reionization process starting at z = 7.1 and going to higher redshifts (to z~12).
"""

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

Dx = 1200. #distance in comoving mpc
N = 1024
print "loading ion"
I9 = n.fromfile('1024_1200Mpc_lc_z7.1.ion',dtype='<f8',count=N**2)
print "loading dens"
M9 = n.fromfile('1024_1200Mpc_lc_z7.1.dens',dtype='<f8',count=N**2)
I9.shape = M9.shape = (N,N)
C = lambda z: 26*n.sqrt((1+z)/10)

X,R = n.mgrid[:N,:N]
R = DM(7.1) + R*Dx/N

Z = z(R)
T = (1+M9)*(1-I9)*C(Z)
X *= Dx/N
#p.figure(facecolor='k',edgecolor='k')
p.figure()
Theta = n.zeros((N,N))
for i in range(N):
    for j in range(N):
        Theta[i,j] = r2theta(X[i,j],Z[i,j]) * 180/n.pi #compute theta in degrees
#ax = p.axes(axisbg='k')
#for line in ax.yaxis.get_ticklines(): 
#    # line is a matplotlib.lines.Line2D instance 
#    line.set_color('w') 
##    line.set_markersize(25) 
##    line.set_markeredgewidth(3) 
#for line in ax.xaxis.get_ticklines(): 
#    line.set_color('w') 
p.imshow(T,cmap='PuBu',extent=(7.1,12,0,8),vmin=0,aspect='auto')
#p.pcolor(Theta,Z,T,cmap='gist_heat',vmin=0)
cb = p.colorbar(orientation='horizontal')
cb.set_label('HI spin temp [mK]')
p.ylabel('degrees')
p.xlabel('z')

p.show()


