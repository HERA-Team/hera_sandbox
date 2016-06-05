#!/usr/bin/env python
#
#  make_sim_figure_v1.py
#  
#
#  Created by Danny Jacobs on 10/1/10.
#  PAPER Project
#
import matplotlib
matplotlib.use('Agg')
from pylab import *
import aipy as a, numpy as n,math as m
import sys, optparse
from cosmo_units import *
import matplotlib as mpl

#make an inverse cosmology model
print "loading cosmology with O_l = %f, O_M = %f, Ho = %f"%(O_l,O_M,Ho)
zi = n.linspace(6,15)
D_i = [DM(z) for z in zi]
z_p = n.polyfit(D_i,zi,4)
z = n.poly1d(z_p)
Dx = 1200.
ds = 6
N = 1024

dtheta = n.degrees(r2theta(Dx/(1024/ds),8))
T_obs = n.fromfile(sys.argv[-1],dtype='<f8')
Xs = n.sqrt(T_obs.size/N)
ds = N/Xs
T_obs.shape = (N/ds,N/ds,N)
D = DM(7.1) + Dx*n.arange(N)/N
print "doing grid cosmology"
Z = n.array([z(d) for d in D])

zmaxs = n.linspace(7.1,10,num=(10-7.1)/0.1)
#compute the image rms for various redshift bin sizes
z_min = 7.
figure()
Trms = []
Tmm = []
for z_max in zmaxs:
    myzs = n.argwhere(n.logical_and(z_min<Z,Z<z_max)).squeeze()
    z_min_i,z_max_i  = (myzs.min(),myzs.max())
    print "integrating %f < z < %f (%d < z_i < %d) "%(z_min,z_max,z_min_i,z_max_i)
    fc = n.mean([f21/(1+z_min),f21/(1+z_max)])
    BW = n.abs(f21/(1+z_min)-f21/(1+z_max))
    print " central frequency =  %6.3f MHz, BW = %3.2f"%(fc/1e6,BW/1e6)
    myT = n.mean(T_obs[:,:,z_min_i:z_max_i],axis=2)
    #myT /= n.mean(myT)
    Trms.append(n.std(myT))
    Tmm.append(myT.max() - myT.min())
#    imshow(myT,cmap=get_cmap('PuBu'),extent=(0,8,0,8),vmin=0,aspect='auto')
#    cb = colorbar(orientation='horizontal')
#    cb.set_label('HI spin temp [mK]')
#    ylabel('deg')
#    xlabel('deg')

subplot(211)
plot(zmaxs,Trms,label='image rms',ls=':',color='k')
plot(zmaxs,Tmm,label='max - min',color='k')
xlabel('redshift')
ylabel('mK')
subplot(212)
BWs = n.abs(f21/(1+zmaxs) - f21/(1+z_min))/1e6
plot(BWs,Trms,ls=':',color='k',label='image rms')
plot(BWs,Tmm,color='k',label='max - min')
xlabel('BW [MHz]')
ylabel('mK')
#filling_factor = 7.**2/((50.*100**2) + 62*(750**2-100**2))#core + ring
filling_factor = (7.**2/(50.*100**2)) #just the core
print "calculating noise using core and outer ring with 112 tiles total"
print "using filling factor = ",filling_factor
hours = 1000
print "at %4.0f hours integration"%hours
Tnoise = 1/filling_factor *  250 / n.sqrt(hours * 3600 * 2)
plot(BWs,Tnoise/n.sqrt(BWs*1e6),'k',lw=3,label='1000 hours with MWA')
legend()
savefig('bw_vs_rms.png')
