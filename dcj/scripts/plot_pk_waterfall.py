#! /usr/bin/env python
"""
Plot a z vs kpl waterfall using the npz files output by pspec_plot_k3pk.py
plot_pk_waterfall.py *npz
"""
from pylab import *
import numpy as n,sys,os

files = sys.argv[1:]
P = []
for file in files:
    F = n.load(file)
    P.append(F['pk'])
kpl = F['kpl']
P = n.array(P)
chans = n.array([int(F.split('/')[1].split('_')[0]) for F in files])
imshow(n.log10(n.abs(P)),aspect='auto',vmin=5,vmax=8.5,extent=(kpl.min(),kpl.max(),chans.max(),chans.min()))
xlabel('$k_\parallel$ [$hMpc^{-1}$]')
ylabel('chan [df=500kHz]')
colorbar()
show()
