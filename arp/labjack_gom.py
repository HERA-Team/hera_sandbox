#! /usr/bin/env python
import numpy as n, aipy as a
import sys

def grid_jd(jds, temps, binsize=120):
    jdbin = binsize * a.ephem.second
    nbins = int((jds[-1] - jds[0]) / jdbin)
    wgts,bins = n.histogram(jds, bins=nbins)
    dats,bins = n.histogram(jds, weights=temps, bins=nbins)
    return dats / wgts, bins

def gain(H, T):
    return 10**(H*(T-T0)/10.)

# GoM E thermistor is tied to a block between the aluminum plates, simulating load temperature.
# Ant58 thermistor is tide to balun casing of Antenna 58.
T_RX, T_ANT58, T_GOM_E = 0,1,2
H_BALUN = -0.024    # dB/K
H_CABLE = -0.018    # dB/K
H_RECVR = -0.045    # dB/K
T0 = 300. # K

lines = [L.split() for f in sys.argv[1:] for L in open(f).readlines()]
jds = n.array([a.phs.ephem2juldate(a.ephem.date(' '.join(L[:2]))) for L in lines])
jds -= 2 * a.ephem.hour # Fix for local time being 2 hrs ahead of UTC
temps = n.array([map(float,L[2:]) for L in lines])
t_rx, bins = grid_jd(jds, temps[:,T_RX])
t_ant58, bins = grid_jd(jds, temps[:,T_ANT58])
t_gom_e, bins = grid_jd(jds, temps[:,T_GOM_E])
jds = 0.5 * (bins[:-1] + bins[1:])

g = gain(H_RECVR, t_rx) * gain(H_BALUN, t_ant58) * gain(H_CABLE, t_ant58)
calgain = n.average(g.compress(n.where(n.abs(jds - 2455747.58) < a.ephem.hour, 1, 0)))
g /= calgain

jds -= 2455747

import pylab as p

p.subplot(211)
p.plot(jds, t_rx, label='recvr')
p.plot(jds, t_ant58, label='ant58')
p.plot(jds, t_gom_e, label='gom_e')
p.legend()

p.subplot(212)
p.plot(jds, g)

p.show()
