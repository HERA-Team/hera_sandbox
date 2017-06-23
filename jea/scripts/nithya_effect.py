# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:56:50 2016

@author: jaguirre
"""

import numpy as np
import healpy as hp
from astropy import units as u
from astropy import constants as c
import pylab as plt

def VisualizeAlm(alm,figno=1,max_l=None):
    """ Visualize a healpy a_lm vector """
    lmax = hp.Alm.getlmax(f_lm.size)
    l,m = hp.Alm.getlm(lmax)
    mag = np.zeros([lmax+1,lmax+1])
    phs = np.zeros([lmax+1,lmax+1])
    mag[m,l] = np.abs(alm)
    phs[m,l] = np.angle(alm)
    cl = hp.alm2cl(alm)
    # Decide the range of l to plot
    if max_l != None:
        max_l = (max_l if (max_l <= lmax) else lmax)
    else:
        max_l = lmax 
    print max_l
    plt.figure(figno)
    plt.clf()
    plt.subplot(211)
    plt.imshow(mag[0:max_l,0:max_l],interpolation='nearest',origin='lower')
    plt.colorbar()
    plt.subplot(212)
    plt.imshow(phs[0:max_l,0:max_l],interpolation='nearest',origin='lower')
    plt.colorbar()
    # plt.subplot(313)
    #plt.semilogy(cl[0:max_l])
    return {'mag':mag,'phs':phs,'cl':cl}
    
#%%

# re-inventing the wheel
def calc_delay(n,df):
    # Convert from MHz to Hz and tau back to ns
    tau = np.fft.fftshift(np.fft.fftfreq(n,d=df*1e6)*1e9)
    return tau

def hp2full_a_lm(a_lm,lmax):
    l,m = hp.Alm.getlm(lmax)
    assert (len(a_lm) == len(l))
    A_lm = {}
    for L in np.arange(lmax+1):
       Mk = {}
       for M in np.arange(-L,L+1):
           i = hp.Alm.getidx(lmax,L,np.abs(M))
           if M < 0:
               a = np.power(-1,M)*np.conj(a_lm[i])
           else:
               a = a_lm[i]
           Mk[M] = a
       A_lm[L] = Mk
    return A_lm

nside = 128
lmax = 3*nside-1
ell = np.arange(0,lmax+1)
l,m = hp.Alm.getlm(lmax)
npix = hp.nside2npix(nside)
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)
below_horizon = np.where(theta > np.pi/2.)[0]

nu = np.linspace(100,200,num=203)*u.MHz
dnu = np.median(np.diff(nu))
tau = calc_delay(len(nu),df=dnu)

## Instrument parameters
lmbda0 = 2.
D = 14.
fwhm = 1.22*lmbda0/D
sigma = fwhm/2.35
b0 = 29.2
tau0 = (b0*u.m/c.c).to(u.ns)

## Construct beam
sidelobe_level = 3e-3
beam = np.exp(-np.power(theta,2)/(2.*np.power(sigma,2)))
beam += sidelobe_level
#beam_nu += sidelobe_level
#beam = np.power(np.sin(theta),2)*np.power(np.sin(phi),16)
beam /= beam.max()
beam_nu = np.outer(beam,np.ones_like(nu))

## Construct fringe
bvec = np.array([0,1,0])*b0
b = np.outer(bvec,np.ones(npix))*u.m
s = np.array(hp.pix2vec(nside,ipix))
bdots = np.sum(b*s,axis=0)
fringe = np.exp(-2.*np.pi*1j*np.outer(bdots.value,nu.to(u.Hz).value)/c.c.value)
if True:
    fringe[below_horizon] = 0.
    beam[below_horizon] = 0.

## V_G(nu) = s_00(u) T(nu) 
# Average compensates for any nside changes, normalized to beam solid angle (?)
omega_nu = 4.*np.pi*beam_nu.sum(axis=0)/npix
T_nu = np.average(beam_nu*fringe,axis=0)
window = np.hanning(len(T_nu))
dtrans_T_nu = np.fft.fftshift(np.fft.fft(window*T_nu))
Ptau_T_nu = np.abs(dtrans_T_nu)

f_lm = hp.map2alm(fringe[:,100])
C_f = hp.alm2cl(f_lm)
a_lm = hp.map2alm(beam_nu[:,100])
C_a = hp.alm2cl(a_lm)

C_fa = hp.alm2cl(f_lm*a_lm)

figno = 5
title = 'Gah'#'Sin^16(phi)'
filename = ''

plt.figure(figno)
plt.clf()
plt.subplot(231)
plt.plot(nu,T_nu.real,'b')
plt.plot(nu,T_nu.real-T_nu.mean(),'b--')
plt.plot(nu,window*T_nu.real,'g')
plt.plot(nu,T_nu.imag,'r')

plt.subplot(232)
plt.plot(nu,T_nu.real/omega_nu)
plt.title(title)

plt.subplot(233)
plt.semilogy(tau,Ptau_T_nu)
plt.semilogy(-tau0.value*np.array([1.,1.]),[1e-2,1e5],'r--')
plt.semilogy(tau0.value*np.array([1.,1.]),[1e-2,1e5],'r--')

plt.subplot(234)
plt.plot(C_f/C_f.max(),label='C_l,fringe')
plt.plot(C_a/C_a.max(),label='C_l,beam')
plt.xlim([0,150])
#plt.ylim([1e-8,1e1])
plt.legend()

plt.subplot(235)
plt.plot(ell,C_fa/C_fa.max())
plt.xlim([0,30])

hp.orthview((beam_nu*fringe)[:,100].real,rot=[0,90],half_sky=True,sub=236)
hp.graticule()


    
    
    
    
    
    
    
    
    