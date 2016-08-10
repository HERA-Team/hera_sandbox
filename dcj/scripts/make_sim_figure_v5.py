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
N = 1024
Dx = 1200.

DZ = 2

#make an inverse cosmology model
print "loading cosmology with O_l = %f, O_M = %f, Ho = %f"%(O_l,O_M,Ho)
zi = n.linspace(6,15)
D_i = [DM(z) for z in zi]
z_p = n.polyfit(D_i,zi,4)
z = n.poly1d(z_p)
print "reading T file"
T = n.fromfile('1024_1200Mpc_lc_z7.1.temp',dtype='<f8')
T.shape = (N,N,N)
N = T.shape[2]
r_i = n.arange(N)
r = DM(7.1) + r_i*Dx/N
zs = z(r)
Z_steps = n.arange(zs.min(),zs.max(),DZ)

p.figure()
U = n.ceil(n.sqrt(len(Z_steps)))
Q = n.ceil(len(Z_steps)/U)
nfigs = U*Q
print "plotting a %d x %d figure"%(U,Q)
for i,Z in enumerate(Z_steps):
    print "plotting z = ",Z
    ax = p.subplot(U,Q,i+1)
    z_slice = n.argwhere(n.logical_and(Z<zs,zs<(Z+DZ))).squeeze()
    plt  = p.imshow(n.mean(T[:,:,z_slice],axis=2),cmap='PuBu',extent=(0,8,0,8),vmin=0,vmax=16,aspect='auto')
    p.text(0.1,0.8,str(Z),transform=ax.transAxes)
#p.pcolor(Theta,Z,T,cmap='gist_heat',vmin=0)
    p.ylabel('degrees')
    p.xlabel('z')
#cb = plt.colorbar(orientation='horizontal')
#cb.set_label('HI spin temp [mK]')

p.show()


