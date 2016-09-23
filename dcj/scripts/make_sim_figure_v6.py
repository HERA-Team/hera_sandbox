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


T_obs = n.fromfile(sys.argv[-1],dtype='<f8')
Xs = n.sqrt(T_obs.size/N)
ds = N/Xs
dtheta = n.degrees(r2theta(Dx/(1024/ds),8))

#T_obs = n.fromfile('1024_1200Mpc_lc_z7.1.temp_smoothsmall',dtype='<f8')

T_obs.shape = (N/ds,N/ds,N)
D = DM(7.1) + Dx*n.arange(N)/N
print "doing grid cosmology"
Z = n.array([z(d) for d in D])



print "plotting"
imgfig = figure(figsize=(16,8))
#start the imaging plot
subplot(231)
imshow(T_obs[:,:,100],cmap=get_cmap('PuBu'),extent=(0,8,0,8),vmin=0,vmax=7,aspect='auto')
#p.pcolor(Theta,Z,T,cmap='gist_heat',vmin=0)
cb = colorbar(orientation='horizontal')
cb.set_label('HI spin temp [mK]')
ylabel('deg')
xlabel('deg')
def plot_zavg_image(z_min,z_max):
    myzs = n.argwhere(n.logical_and(z_min<Z,Z<z_max)).squeeze()
    z_min_i,z_max_i  = (myzs.min(),myzs.max())
    print "integrating %f < z < %f (%d < z_i < %d) "%(z_min,z_max,z_min_i,z_max_i)
    fc = n.mean([f21/(1+z_min),f21/(1+z_max)])
    BW = n.abs(f21/(1+z_min)-f21/(1+z_max))
    print " central frequency =  %6.3f MHz, BW = %3.2f"%(fc/1e6,BW/1e6)
    myT = n.mean(T_obs[:,:,z_min_i:z_max_i],axis=2)
    #myT /= n.mean(myT)
    imshow(myT,cmap=get_cmap('PuBu'),extent=(0,8,0,8),aspect='auto')
    cb = colorbar(orientation='horizontal')
    cb.set_label('HI spin temp [mK]')
    ylabel('deg')
    xlabel('deg')


#start the correlation width figure
corrfig = figure(figsize=(16,8))
x_slice = 6.0
print "selecting slice at x=%3.1f deg"%x_slice
thetas = n.arange(0,dtheta*N/ds,dtheta)
myx = n.int(n.median(n.argwhere(n.logical_and((x_slice-dtheta*2)<thetas,thetas<(x_slice+dtheta*2))).squeeze()))
print "index = %d"%myx

subplot(231)
imshow(T_obs[:,myx,:],cmap=get_cmap('PuBu'),extent=(7.1,12,0,8),vmin=0, vmax=10,aspect='auto')
#pcolor(Theta,Z,T,cmap='gist_heat',vmin=0)
cb = colorbar(orientation='horizontal')
cb.set_label('HI spin temp [mK]')
ylabel('degrees')
xlabel('z')



def plot_width_image(z_min,z_max):
    myzs = n.argwhere(n.logical_and(z_min<Z,Z<z_max)).squeeze()
    z_min_i,z_max_i  = (myzs.min(),myzs.max())
    print "integrating %f < z < %f (%d < z_i < %d) "%(z_min,z_max,z_min_i,z_max_i)
    fc = n.mean([f21/(1+z_min),f21/(1+z_max)])
    BW = n.abs(f21/(1+z_min)-f21/(1+z_max))
    print " central frequency =  %6.3f MHz, BW = %3.2f"%(fc/1e6,BW/1e6)
    #compute the redshift direction correlation
    myT = T_obs[:,:,z_min_i:z_max_i]/n.mean(T_obs[:,:,z_min_i:z_max_i])
    FT = n.fft.fft(myT,axis=2)
    P2 = FT*n.conj(FT)
    nz = z_max_i - z_min_i # # of redshift samples
    dz = n.abs(z_max - z_min)/float(nz) #redshift unit
    Dz = n.abs(z_max - z_min)
    A = n.tile(n.fft.fftfreq(nz,d=dz),(N/ds,N/ds,1)) #inverse redshift distance axis
    
    #compute and plot an image of the width
    W = n.sqrt(n.abs(n.sum(A**2*P2,axis=2)/n.sum(P2,axis=2)))
    print W.min(),n.median(W),W.max()
    imshow(W,cmap=get_cmap('PuBu'),extent=(0,8,0,8),aspect='auto',vmin=0.9,vmax=2)
    
    #fix up the plot
    A_ticks = n.linspace(1/2.,2,num=4) * n.median(W) #inverse redshift ticklabels
    A_ticks = n.linspace(0.9,2,num=4)
    Zi_ticks = 1/A_ticks #redshift correlation length ticklabels
    cb = colorbar(orientation='horizontal',ticks=A_ticks)
    cb.set_label('correlation length [dz]')
    cb.set_ticklabels(map(lambda x: "%3.2f"%x,Zi_ticks))
    ylabel('deg')
    xlabel('deg')


myzs = [[7.0,8.0],
#        [7.25,7.75],
        [8.0,9.0],
#        [7.75,8.25],
        [9.0,10.0],
        [10.0,11.],
        [10.5,11.5]]
for i,zs in enumerate(myzs):
    figure(imgfig.number)
    ax = subplot(2,3,i+2)
    plot_zavg_image(zs[0],zs[1])
    text(0.6,0.8,'$%3.1f<z<%3.1f$'%(zs[0],zs[1]),transform=ax.transAxes,color='y')

    figure(corrfig.number)
    ax = subplot(2,3,i+2)
    plot_width_image(zs[0],zs[1])
    text(0.6,0.8,'$%3.1f<z<%3.1f$'%(zs[0],zs[1]),transform=ax.transAxes,color='k')
    
savefig('smooth_temps_correlation_width.png')

figure(imgfig.number)
savefig('smooth_temps_experiment.png')


#figure(figsize=(16,4))
#x_slice = 5.0
#print "selecting slice at x=%3.1f deg"%x_slice
#thetas = n.arange(0,dtheta*N/ds,dtheta)
#myx = n.int(n.median(n.argwhere(n.logical_and((x_slice-dtheta*4)<thetas,thetas<(x_slice+dtheta*4))).squeeze()))
#print "index = %d"%myx
#subplot(131)
#imshow(T_obs[:,myx,:],cmap=get_cmap('PuBu'),extent=(7.1,12,0,8),vmin=0, vmax=5,aspect='auto')
##pcolor(Theta,Z,T,cmap='gist_heat',vmin=0)
#cb = colorbar(orientation='horizontal')
#cb.set_label('HI spin temp [mK]')
#ylabel('degrees')
#xlabel('z')
#
#subplot(132)
#z_min,z_max = (7,8)
#myzs = n.argwhere(n.logical_and(z_min<Z,Z<z_max)).squeeze()
#z_min_i,z_max_i  = (myzs.min(),myzs.max())
#print "integrating %f < z < %f (%d < z_i < %d) "%(z_min,z_max,z_min_i,z_max_i)
#fc = n.mean([f21/(1+z_min),f21/(1+z_max)])
#BW = n.abs(f21/(1+z_min)-f21/(1+z_max))
#print " central frequency =  %6.3f MHz, BW = %3.2f"%(fc/1e6,BW/1e6)
#FT = n.fft.fft(T_obs[:,myx,z_min_i:z_max_i],axis=1)
#print n.abs(FT).min(), n.abs(FT).max()
#imshow(n.log(n.abs(FT*n.conj(FT))),cmap=get_cmap('PuBu'),extent=(z_min,z_max,0,8),aspect='auto')
#cb = colorbar(orientation='horizontal')
#cb.set_label('HI spin temp [mK]')
#ylabel('deg')
#xlabel('deg')
#
#
#subplot(133)
#z_min,z_max = (7.5,8.5)
#myzs = n.argwhere(n.logical_and(z_min<Z,Z<z_max)).squeeze()
#z_min_i,z_max_i  = (myzs.min(),myzs.max())
#print "integrating %f < z < %f (%d < z_i < %d) "%(z_min,z_max,z_min_i,z_max_i)
#FT = n.fft.fft(T_obs[:,myx,z_min_i:z_max_i],axis=1)
#imshow(n.log(n.abs(FT*n.conj(FT))),cmap=get_cmap('PuBu'),extent=(z_min,z_max,0,8),aspect='auto')
##semilogy(n.abs(FT*n.conj(FT))[0,:])
##pcolor(Theta,Z,T,cmap='gist_heat',vmin=0)
#cb = colorbar(orientation='horizontal')
#cb.set_label('HI spin temp [mK]')
#ylabel('degrees')
#xlabel('z')

#savefig('smooth_temps_correlation.png')



