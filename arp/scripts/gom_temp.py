#!/usr/bin/env python
import aipy as a, numpy as n, sys

def to_dB(dat): return 10*n.log10(dat)
def from_dB(dat): return 10**(dat / 10)

dat = []
gom = []
for arg in sys.argv[1:]:
    print 'Reading', arg
    uv = a.miriad.UV(arg)
    uv.select('antennae', 1, 1)
    for p,d,f in uv.all(raw=True):
        (crd, t, (i,j)) = p
        dat.append([t, uv['t_recvr'], uv['t_cable'], uv['t_balun'], uv['t_load']])
        gom.append(d)

dat = n.array(dat)
gom = n.median(n.array(gom).real, axis=1)
gom_c = gom / (dat[:,4]/300)  # calibrate to load temp of 300K
gom_c_dB = to_dB(gom_c)

def fit_1H(H, T):
    return n.std(gom_c_dB - H[0] * T)

def fit_2H(Hs, T_c, T_b):
    H_c, H_b = Hs
    return n.std(gom_c_dB - H_c * T_c - H_b * T_b)

H_c1 = a.optimize.fmin(fit_1H, [0], args=(dat[:,2],))
H_b1 = a.optimize.fmin(fit_1H, [0], args=(dat[:,3],))
H_c,H_b = a.optimize.fmin(fit_2H, [0,0], args=(dat[:,2],dat[:,3]))

print 'H_c1:', H_c1
print 'H_b1:', H_b1
print 'H_c/H_b:', H_c, H_b

import pylab as p

p.subplot(211)
p.plot(dat[:,0], dat[:,1], '-', label='recvr')
p.plot(dat[:,0], dat[:,2], '-', label='cable')
p.plot(dat[:,0], dat[:,3], '-', label='balun')
p.plot(dat[:,0], dat[:,4], '-', label='load')
p.legend()

p.subplot(212)
p.plot(dat[:,0], gom, '-', label='gom')
p.plot(dat[:,0], gom_c, '-', label='gom/load')
p.plot(dat[:,0], gom_c / from_dB(H_c1 * (dat[:,2] - n.average(dat[:,2]))), 
    '-', label='cable')
p.plot(dat[:,0], gom_c / from_dB(H_b1 * (dat[:,3] - n.average(dat[:,3]))), 
    '-', label='balun')
p.plot(dat[:,0], gom_c / from_dB(H_c * (dat[:,2] - n.average(dat[:,2])) + H_b * (dat[:,3] - n.average(dat[:,3]))), 
    '-', label='both')
#p.legend()
p.show()
