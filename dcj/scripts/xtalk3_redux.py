#!/usr/bin/env python
#
#  xtalk3_redux.py
#  
#
#  Created by Danny Jacobs on 8/4/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pickle
from pylab import *
from mpl_toolkits.axes_grid import make_axes_locatable
from smooth import smooth
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--snap', dest='snap', action='store_true',
    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

"""
Loads .xtalk files and generates a smooth, time dependent model to subtract from data
"""
I = n.complex(0,1)
data ={}
for file in args:
    t = float(file[:-6])
    print "loading ",file
    data[t] = pickle.load(open(file,'r'))
ts = n.sort(data.keys())
#reorder by bl then time
bls = data[data.keys()[0]].keys()
xtalk = {}
for bl in bls:
    D = n.array([data[t][bl] for t in ts])
    xtalk[bl] = n.ma.masked_where(D==-n.Inf,n.ma.masked_equal(D,0))
pickle.dump(xtalk,open('xtalk.pkl','w'))
model = {}
t = ts
f = n.arange(xtalk[bl].shape[1])



#ax1  = subplot(211)
#ax2 = subplot(212)
#for bl in xtalk:
bl=1547

XT_t = n.average(xtalk[bl],axis=1)
XT_f = n.average(xtalk[bl],axis=0)
chans = n.argwhere(XT_f!=0)
f = f[chans].squeeze()
print f.shape
F,T = n.meshgrid(f,t)
XT_f = XT_f[chans].squeeze()

n.save('tdump',n.average(xtalk[bl],axis=1))
n.savez('fdump',XT_f=XT_f,f=f)

print "fitting baseline %d_%d"%a.miriad.bl2ij(bl)
#t_amp = n.poly1d(n.polyfit(t,n.abs(XT_t),2))
#t_ang = n.poly1d(n.polyfit(t,n.angle(XT_t),2))
#f_amp = n.poly1d(n.polyfit(f,n.abs(XT_f),2))
#f_ang = n.poly1d(n.polyfit(f,n.angle(XT_f),2))
#t_model = t_amp(t)*n.exp(I*t_ang(t))
#f_model = f_amp(f)*n.exp(I*f_ang(f))

N = 10 #number of free parameters in smooth
flag_err = 0.5 #fractional error above which to flag
t_window_len = N-len(XT_t)
f_window_len = N-len(XT_f)
t_model = smooth(n.abs(XT_t),window_len=t_window_len)*\
    n.exp(I*smooth(n.angle(XT_t),window_len=t_window_len))
t_model = smooth(n.abs(XT_f),window_len=f_window_len)*\
    n.exp(I*smooth(n.angle(XT_f),window_len=f_window_len))
model[bl] = n.outer(t_model,f_model)/n.average(xtalk[bl])
Model_Angle = n.angle(model[bl])
Model_Amp = n.abs(model[bl])

Xtalk_Angle = n.angle(xtalk[bl])
flags = n.argwhere(n.logical_or(model[bl])

pickle.dump(model[bl],open('test.pkl','w'))

bl=1547

figure(10)
subplot(311)
pcolor(F,T,n.log10(n.abs((xtalk[bl][:,f]-model[bl])/xtalk[bl][:,f])))
colorbar()

ax1 = subplot(312)
ax1.plot(t,n.abs(XT_t),'.')
ax1.plot(t,t_amp(t))
xlabel('time [int]')
ylabel('amplitude')
ax2 = subplot(313)

ax2.semilogy(f,n.abs(XT_f),'.')
ax2.semilogy(f,f_amp(f))
xlabel('frequency [chan]')
ylabel('amplitude')

figure(11)
subplot(311)
pcolor(F,T,n.angle((xtalk[bl][:,f]-model[bl])/xtalk[bl][:,f]))
colorbar()

ax1 = subplot(312)
ax1.plot(t,n.angle(XT_t),'.')
ax1.plot(t,t_ang(t))
xlabel('time [int]')
ylabel('phase [r]')
ax2 = subplot(313)

ax2.plot(f,n.angle(XT_f),'.')
ax2.plot(f,f_ang(f))
xlabel('frequency [chan]')
ylabel('phase [r]')

show()

    