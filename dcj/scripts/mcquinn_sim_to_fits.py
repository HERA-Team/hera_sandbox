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
from astropy.io import fits
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

o= optparse.OptionParser()
o.add_option('-z',type='float',default=8.,help='center redshift [default=8]')
o.add_option('--dz',default=0.5,help='redshift window [default=0.5]')
#a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])
zmax = opts.z+opts.dz/2
zmin = opts.z-opts.dz/2
#make an inverse cosmology model
print "loading cosmology with O_l = %f, O_M = %f, Ho = %f"%(O_l,O_M,Ho)
zi = n.linspace(6,15)
D_i = [DM(z) for z in zi]
z_p = n.polyfit(D_i,zi,4)
z = n.poly1d(z_p)

Dx = 1200. #distance in comoving mpc
N = 1024
print "creating coordinate system"
#X,Y,R = n.mgrid[:N,:N,:N]
#del(X);del(Y)
#R = DM(7.1) + R*Dx/N
r = DM(7.1) + n.arange(N)*Dx/N
zs = z(r)


Nz = n.logical_and(zs>zmin,z<zmax).sum()
#Z = z(R);del(R)
J = n.argwhere(n.logical_and(zs>zmin,zs<zmax))
zs = zs[J].squeeze()

print "loading ion"
I9 = n.fromfile('1024_1200Mpc_lc_z7.1.ion',dtype='<f8')
I9.shape = (N,N,N)
I9=I9[:,:,J]
print "Average ionization fraction at z=%3.1f : %3.2f"%(opts.z,n.mean(I9)) 
print "loading dens"
M9 = n.fromfile('1024_1200Mpc_lc_z7.1.dens',dtype='<f8')
M9.shape = (N,N,N)
M9 = M9[:,:,J]

#Z = Z[J]
C = lambda z: 26*n.sqrt((1+z)/10)
print "done loading`"
print "calculatings"
print zs
print f21/(1+zs.max()-n.diff(zs)[-1])-f21/(1+zs.max())
Cs = C(zs).squeeze()
T = (1+M9)
del(M9)
T *= (1-I9)
del(I9)
print Cs.shape,T.shape
T = T.squeeze()
T *= Cs
print "writing file"
RA=0
DEC=-30
RES=Dx/N/DM(opts.z)*180/n.pi
T = T.T
hdu = fits.PrimaryHDU(T)
print T.shape
hdu.header['CTYPE1']='RA---SIN'
hdu.header['CRPIX1']=T.shape[2]/2 #RApix
hdu.header['CRVAL1']=0.0 #center RA
hdu.header['CDELT1']=RES
hdu.header['CTYPE2']='DEC--SIN'
hdu.header['CRPIX2']=T.shape[1]/2 #DECpix
hdu.header['CRVAL2']=-30.0 #DEC
hdu.header['CDELT2']=RES
hdu.header['CTYPE3']='FREQ'
hdu.header['CRPIX3']=0            #FREQ
hdu.header['CRVAL3']=f21/(1+zs.max())
hdu.header['CDELT3']=f21/(1+zs.max()-n.diff(zs)[-1])-f21/(1+zs.max())
hdu.writeto('1024_1200Mpc_lc_z%3.1f.fits'%opts.z,clobber=True)
#T.tofile('1024_1200Mpc_lc_z7.1.temp')
sys.exit()
#X *= Dx/N
#p.figure(facecolor='k',edgecolor='k')
p.figure()
#Theta = n.zeros((N,N))
#for i in range(N):
#    for j in range(N):
#        Theta[i,j] = r2theta(X[i,j],Z[i,j]) * 180/n.pi #compute theta in degrees
#ax = p.axes(axisbg='k')
#for line in ax.yaxis.get_ticklines(): 
#    # line is a matplotlib.lines.Line2D instance 
#    line.set_color('w') 
##    line.set_markersize(25) 
##    line.set_markeredgewidth(3) 
#for line in ax.xaxis.get_ticklines(): 
#    line.set_color('w') 
p.imshow(T[:,:,400],cmap='PuBu',extent=(7.1,12,0,8),vmin=0,aspect='auto')
#p.pcolor(Theta,Z,T,cmap='gist_heat',vmin=0)
cb = p.colorbar(orientation='horizontal')
cb.set_label('HI spin temp [mK]')
p.ylabel('degrees')
p.xlabel('z')

p.show()


