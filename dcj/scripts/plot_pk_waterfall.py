#! /usr/bin/env python
"""
Plot a z vs kpl waterfall using the npz files output by pspec_plot_k3pk.py
plot_pk_waterfall.py *npz
"""
from pylab import *
import numpy as n,sys,os
from capo.cosmo_units import *

files = sort(sys.argv[1:])
P = []
Perr = []
for file in files:
    print file
    F = n.load(file)
    P.append(F['pk'])
    Perr.append(F['err'])
kpl = F['kpl']
chans = n.array([int(F.split('/')[1].split('_')[0]) for F in files])
I = n.argsort(chans)
chans = chans[I]
P = n.array(P)[I]
Perr = n.array(Perr)[I]

freqs = chans/2. + 100
redshifts = n.array(list(set(n.round(f212z(freqs*1e6)))))
zfreqs = f21/(1+redshifts)/1e6
imshow(n.log10(n.abs(P)),aspect='auto',vmin=5,vmax=9,extent=(kpl.min(),kpl.max(),chans.max(),chans.min()))
xlabel('$k_\parallel$ [$hMpc^{-1}$]')
ylabel('chan [df=500kHz]')
colorbar()

#plot a single k// bin
mykpl = 0.35
ki = []
kineg = []
for file in files:
    F = n.load(file)
    ki.append(n.abs(F['kpl'] - mykpl).argmin())
    kineg.append(n.abs(F['kpl'] + mykpl).argmin())
print ki
Pk_at_mykpl =      n.array( [P[i,ki[i]] for i in range(len(ki))])
Pk_err_kpl =       n.array(  [Perr[i,ki[i]] for i in range(len(ki))])
Pk_at_mykpl_neg =  n.array(  [P[i,kineg[i]] for i in range(len(kineg))])
Pk_err_kpl_neg =   n.array(  [Perr[i,kineg[i]] for i in range(len(kineg))])

figure()
ax1 = subplot(111)
#ax1.errorbar(freqs,Pk_at_mykpl,yerr=Pk_err_kpl)
#ax1.errorbar(freqs,Pk_at_mykpl_neg,yerr=Pk_err_kpl_neg)
Pkk = (Pk_at_mykpl+Pk_at_mykpl_neg)/2
Pkk_err = n.sqrt((Pk_err_kpl**2 + Pk_err_kpl_neg**2)/2)
ax1.set_yscale('log',nonposy='clip')
ax1.errorbar(freqs,Pkk,yerr=Pkk_err,fmt='xk')
ax1.set_ylim([10**6,10**9])
xlims = ax1.get_xlim()
ax2 = ax1.twiny()
ax2.set_xlim(ax1.get_xlim())
print zfreqs,redshifts
ax2.set_xticks(zfreqs)
ax2.set_xticklabels(redshifts)

show()
