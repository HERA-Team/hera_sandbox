#!/usr/bin/env python
#
#  RRL_sim.py
#  
#
#  Created by Danny Jacobs on 8/24/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
from pylab import *

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--snap', dest='snap', action='store_true',
    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

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
        if n.isnan(Y_out[i]):
            print "WARNING: empty grid locations.  RM synthesis will suffer! Decrease number of RM channels to avoid this problem."
    return Y_out,L2,wgts
    
N = n.arange(300,400)
Rnu = 3.289841015e15
RRL = Rnu*(1./N**2-1./(N+1.)**2)
c = 3e8
F = n.linspace(100,200,num=4096)
M = n.zeros_like(F)
RRL = RRL[20:]/1e6
rand_RRL = n.random.uniform(100,200,size=len(N))
rand_M = n.zeros_like(F)
for f in RRL:
    C = n.argwhere(n.abs(f-F)==n.min(n.abs(f-F))).squeeze()
    M[C] = 1
for f in rand_RRL:
    C = n.argwhere(n.abs(f-F)==n.min(n.abs(f-F))).squeeze()
    rand_M[C] = 1
figure(2)   
clf() 
for i in range(200,600,50):
    print i
    loglog(n.abs(n.fft.ifft(M[i:i+50])))
figure(3)
clf()
Ncosm = 50
RMtest = 10
subplot(211)
plot(c/(F*1e6),n.ma.masked_where(M==0,M),'.')
plot(c/(F*1e6),rand_M,alpha=0.4)
plot(c/(F*1e6),sin((c/(F*1e6))**2*2*n.pi*RMtest),'0.5')
xlabel('wavelength [m]')

subplot(212)
ML2,L2,W = grid_L2(M,c/(F*1e6),num=len(F)/4)

rand_ML2,L2,rand_W = grid_L2(rand_M,c/(F*1e6),num=len(F)/4)
spike_ML2,L2,spike_W =  grid_L2(sin((c/(F*1e6))**2*2*n.pi*RMtest),c/(F*1e6),num=len(F)/4)

RM = n.fft.fftfreq(len(ML2),d=n.diff(L2)[0])[len(ML2)/2:]
ML2 /= W
rand_ML2 /= rand_W
spike_ML2 /= spike_W
ML2[n.argwhere(W==0)] = 0
rand_ML2[n.argwhere(rand_W==0)] = 0
spike_ML2[n.argwhere(rand_W==0)] = 0

RMsyn = n.fft.ifft(ML2)[len(ML2)/2:]
rand_RMsyn = n.fft.ifft(rand_ML2)[len(ML2)/2:]
spike_RMsyn = n.fft.ifft(spike_ML2)[len(ML2)/2:]

h = plot(RM,n.abs(RMsyn))[0]
plot(RM,n.abs(rand_RMsyn),alpha=0.4)
plot(RM,n.abs(spike_RMsyn),'k')
xlabel('rotation measure [$m^{-2}$]')
draw()