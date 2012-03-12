#!/usr/bin/env python
#
#  angle_sim1.py
#  
#
#  Created by Danny Jacobs on 11/2/09.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,ephem
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

#look in Taylor 09 for huge RM catalog
file = open("/Users/danny/Work/radio_astronomy/stokes/RMCatalogue.txt")
dtext = [('ra','<f8'),('dec','f'),('flux','f'),
    ('pol_flux','f'),('pol_perc','f'),('rm','f')]
def recenter_FFT(F):
    return n.concatenate([F[len(F)/2+1:],F[:len(F)/2+1]])
def inner_nyquist(F):
    return n.concatenate([F[len(F)/4:len(F)*1/2],F[len(F)*1/2:len(F)*3/4]])
def slot(B,v):
    """
    returns the index (j) into B such that B_j>v>B_j+1 or B_j<v<B_j+1
    Assumes values of B define edges of bins, as in those returned by
    numpy.histogram. If v lies outside of B, return False.
    """
    if v>n.max(B) or v<n.min(B): return False
    for i,b in enumerate(B):
        if i==len(B)-1: continue
        elif v>b and v<B[i+1]: return i
        elif v<b and v>B[i+1]: return i 
def grid_L2(Y,L,num=None):
    """
    Transform: Y(L) -> Y'(L**2)
    where num=len(L**2) desired (default num = len(L))
    return: Y',L2,weights    
    """
    if num is None: num=len(L)
    L2 = n.linspace(n.max(L**2),n.min(L**2),num=num) #frequency sorted!
    Y_out = n.zeros_like(L2)
    wgts = n.zeros_like(L2)
    for i in range(1,len(L2)):
        Y_out[i] = n.average(Y[n.where(
            [a and b for (a,b) in zip(L**2<L2[i-1],L**2>L2[i])])[0]])
        wgts[i] = n.sum([a and b for (a,b) in zip(L**2<L2[i-1],L**2>L2[i])])
    return Y_out,L2,wgts
        

taylor = []
for line in file.readlines():
    ra = ephem.hours(':'.join(line[1:12].split()))
    dec = ephem.degrees(':'.join(line[18:30].split()))
    flux = float(line[69:75])
    pol_flux = float(line[91:95])
    pol_perc = float(line[110:115])
    rm = float(line[127:133])
    taylor.append((ra,dec,flux,pol_flux,pol_perc,rm))
taylor = n.array(taylor,dtype=dtext)
file.close()

al = 1.2
#dP = 0.1 #1sigma confusion limit for PAPER [Jy] (single beam flux due to unresolved sources)
dP = 0.08 #7 days, 4 hours/day * 100kHz channels 
num = 2048
f = n.linspace(100e6,200e6,num)
c = 3e8
L = c/f
L2 = n.linspace(n.max(L**2),n.min(L**2),num=len(L))
M = []


print "sources meeting the taylor criteria in the paper band"
print "assuming confusion noise and an average of 5% polarization"
print n.sum(n.where(taylor['pol_flux']/1000*(3/0.2)**0.7>(1*dP),1,0))



for i,src in enumerate(taylor):
    if src['pol_flux']/1000*(1420./150)**0.7>(1*dP):
        print float(i)/len(taylor)
        P = src['pol_flux']/1000*(1/0.21)**(al)
        t = src['rm']*L**2
        U = n.sin(2*t)*P + n.random.normal(0,dP,num)
        Q = n.cos(2*t)*P + n.random.normal(0,dP,num)
#        ta = src['rm']*L2
#        Ua = n.sin(2*ta)*P# + n.random.normal(0,dP,num)
#        Qa = n.cos(2*ta)*P# + n.random.normal(0,dP,num)


        Ug,L2,Uwgts = grid_L2(U,L,num=len(L)/8)
        Qg,L2,Qwgts = grid_L2(Q,L,num=len(L)/8)
        RM_max = len(L2)/(n.max(L2)-n.min(L2))
        RM = n.linspace(RM_max/4,-RM_max/4,num=len(L2)/2)*n.pi
#        S = [n.complex(q,u) for (q,u) in n.vstack([Q,U]).transpose()]
        Sg =[n.complex(q,u) for (q,u) in n.vstack([Qg,Ug]).transpose()] 
#        FTS = n.fft.fft(S)
        FTSg = n.fft.fft(Sg)
#        PS = n.abs(inner_nyquist(recenter_FFT(FTS)))
        PSg = n.abs(inner_nyquist(recenter_FFT(FTSg))) 
        dRM  = (n.max(RM[n.argwhere(PSg>0.1*n.max(PSg))])-\
            n.min(RM[n.argwhere(PSg>0.1*n.max(PSg))]))
        if dRM !=0:
            M.append(
                (i,src['rm'],RM[n.argwhere(PSg\
                ==n.max(PSg))],
                dRM)
                )
    

sys.exit()
#in NVSS, select bright polarized fluxes, sort by pol flux, look for low dec
#using 000017-341028
num=2048
P  = 11.55#mJy polarized flux at 1.4GHz
#P = #test fluxmJy
#dP = 1.73 #NVSS number
dP = P*0.3
L21 = 0.21 #observation wavelength in meters
#L = n.linspace(3,1.5,num)


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