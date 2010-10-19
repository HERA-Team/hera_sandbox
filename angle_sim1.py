#!/usr/bin/env python
#
#  angle_sim1.py
#  
#
#  Created by Danny Jacobs on 11/2/09.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
from scipy.optimize import fmin
#o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
#opts, args = o.parse_args(sys.argv[1:])

def unwrap(chi,f):
    "Unwrap pi/2 phase changes if phase jumps by more than a fraction (f) of pi/2"
    wrap_delta=0
    chi_uw = n.zeros_like(chi)
    for i in range(1,len(chi)):
        if n.abs(chi[i]-chi[i-1])>(n.pi/2*f):wrap_delta+=-chi[i]+chi[i-1]
        chi_uw[i] = chi[i]+wrap_delta 
    return chi_uw


#in NVSS, select bright polarized fluxes, sort by pol flux, look for low dec
#using 000017-341028
num=2048
P  = 11.55#mJy polarized flux at 1.4GHz
#P = #test fluxmJy
#dP = 1.73 #NVSS number
dP = P*0.3
L21 = 0.21
#L = n.linspace(3,1.5,num)
f = n.linspace(100e6,200e6,num)
c = 3e8
L = c/f
al = 1.2
#make a figure showing effect of channel width on polarization wrapping
p.figure(1)
p.clf()
min_chan = 3
p.semilogy(c/L[1:],n.abs(1/(n.diff(30*(L**2))/(n.pi/2)))/min_chan,label='RM=30')
p.semilogy(c/L[1:],n.abs(1/(n.diff(100*(L**2))/(n.pi/2)))/min_chan,label='RM=100')
p.semilogy(c/L[1:],n.abs(1/(n.diff(300*(L**2))/(n.pi/2)))/min_chan,label='RM=300')
p.legend(loc='upper right')
p.xlabel('frequency [Hz]')
p.ylabel('RM SNR')


p.figure(2)
p.clf()
for i,RM in enumerate([30,100,300]):#RM Units r/m^2
    t = RM*L**2
    p.plot(f[1:],n.diff(t)/n.diff(L**2),'k')
    print "RM ang:",RM, n.max(t)
    t -= n.max(t)
    #U = n.random.normal(P,dP)(L/L21)**al
    Q = n.cos(2*t)*n.random.normal(P,dP,num)*(L/L21)**al  #scale to the PAPER band
    U_pristine = n.clip(n.tan(2*t),-5,5)*Q
#    U =n.clip(n.tan(2*t),-5,5)*Q+n.random.normal(0,dP,num) #Clip to avoid tan divergence.
    U = n.sin(2*t)*n.random.normal(P,dP,num)*(L/L21)**al
    chi = 1/2.*n.arctan(U/Q)
    chi_ma = n.ma.masked_array(chi[2:],n.abs(n.diff(n.diff(chi)/n.diff(L**2)))>100)
    p.plot(f[3:],n.diff(chi_ma)/n.diff(L[2:]**2),label=str(RM))
    RM_m = n.ma.average(n.diff(chi_ma)/n.diff(L[2:]**2))
    RM_s = n.ma.sqrt(n.ma.var(chi_ma))
    print "Simulated source:\n P, dP, RM",P,dP,RM
    print "Recovered RM, dRM:",RM_m,RM_s
p.legend(loc='upper left')    
p.xlabel('freq [Hz]')    
p.ylabel('Rotation Measure [r/m^2]')
#chi_uw = unwrap(chi,0.6)
#def angle_min(V):
#    RM = V[0]
#    ang = V[1]
#    print RM
#    return n.abs(n.average(chi_uw-RM*L**2-ang))
#RM,ang = fmin(angle_min,[100,0])
#print "opt RM ang:",RM,ang    

RM = 30
t = RM*L**2
p.plot(f[1:],n.diff(t)/n.diff(L**2),'k')
#print "RM ang:",RM, n.max(t)
t -= n.max(t)
#U = n.random.normal(P,dP)(L/L21)**al
Q = n.cos(2*t)*n.random.normal(P,dP,num)*(L/L21)**al  #scale to the PAPER band
U_pristine = n.clip(n.tan(2*t),-5,5)*Q
U = n.sin(2*t)*n.random.normal(P,dP,num)*(L/L21)**al
chi = 1/2.*n.arctan(U/Q)
chi_ma = n.ma.masked_array(chi[2:],n.abs(n.diff(n.diff(chi)/n.diff(L**2)))>100)
p.figure(3)
p.clf()
p.subplot(211)
p.plot(f[1:],n.diff(chi)/n.diff(L**2))
p.plot(f[1:],n.diff(t)/n.diff(L**2),'k')
p.title('raw differential RM')
p.xlabel('freq [Hz]')
p.subplot(212)
p.plot(f[3:],n.diff(chi_ma)/n.diff(L[2:]**2))
p.plot(f[1:],n.diff(t)/n.diff(L**2),'k')
p.title('cleaned differential RM')
p.xlabel('freq [Hz]')
p.ylabel('RM [r/m^2]')
#
#RM_low = 30 #r/m^2
#t = RM_low*L**2
##U = n.random.normal(P,dP)(L/L21)**al
#Q = n.random.normal(P,dP,num)*(L/L21)**al
#U = n.tan(2*t)*Q+n.random.normal(0,dP,num)
#p.subplot(212)
#p.plot(1/2.*n.arctan(U/Q)/L**2,'.')




p.show()