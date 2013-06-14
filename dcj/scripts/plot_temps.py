#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse,os
from pylab import *
o = optparse.OptionParser()
o.add_option('--H_balun',dest='H_balun', type='float', default=-0.024,
    help='Temperature coefficient (dB/K) of balun gain.  Default -0.024')
o.add_option('--H_cable',dest='H_cable', type='float', default=-0.018,
    help='Temperature coefficient (dB/K) of cable gain.  Default -0.018')
opts,args = o.parse_args(sys.argv[1:])


def gain(H, T): return 10**(H*(T - 300)/10)

def dB(x): return 10*n.log10(x)

def mfunc(uv, p, d, f):
    crd,t,(i,j) = p
    g = gain(opts.H_balun, uv['t_balun']) * gain(opts.H_cable, uv['t_cable'])
    if i in goms: g *= n.sqrt(uv['t_load'] / opts.T0)
    if j in goms: g *= n.sqrt(uv['t_load'] / opts.T0)
    return p, d/g, f
t_load = []
t_balun = []
t_cable = []
times = []
lst = []
for filename in args:
    uvi = a.miriad.UV(filename)
    print filename
    sys.stdout.flush()
    for (uvw,t,(i,j)),d in uvi.all():
        t_load.append(uvi['t_load'])
        t_balun.append(uvi['t_balun'])
        t_cable.append(uvi['t_cable'])
        times.append(t)
        lst.append(uvi['lst'])
figure()
t_load = n.array(t_load)
t_balun = n.array(t_balun)
t_cable = n.array(t_cable)
times = n.array(times)
plot(times-times.min(),t_load,'-',label='Load',color='0.5')
plot(times-times.min(),t_balun,'-k',label='Balun',lw=5)
plot(times-times.min(),t_cable,'-b',label='Cable')
xlabel('time since %f'%times.min())
ylabel('temp [K]')
twinx()
plot(times-times.min(),dB(gain(opts.H_balun,t_balun)),':k')
plot(times-times.min(),dB(gain(opts.H_cable,t_cable)),':b')
ylabel('gain [dB]')
show()
