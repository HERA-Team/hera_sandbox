#! /usr/bin/env python
"""
Plot the output of uv_rms.py (uv_rms.npz) 

Takes multiple npz files and uses the label variable in the npz to detirmine line style

"""

from pylab import *
import numpy as n
from glob import glob
import matplotlib
from capo.cosmo_units import *
TMAX=1e2 #maximum plot in mK
TMIN=1#minimum plot in mK
w_earth = 2*pi*1.16e-5 #sidereal rotation rate of earth in rad/s
#matplotlib.rcParams.update({'font.size':16})

FILENAME_LABELS=False #label by filenames instead of by npz label parameter
#LOAD AND PLOT THE output from uv_rms.py
files = sys.argv[1:]
linestyles = {'dbl':'b','dbl,fr':'b','df':'m','df,fr':'m'}
ax1 = subplot(111)
ax1.set_yscale('log')
for file in files:
    F = n.load(file)
    freqs = F['freqs']
    flux_scale = (363.4/424) * (freqs/0.160)**(-0.76 + .95) #copied from noise_vs_fq.py by Parsons
    #sets the psa32 data set to the Jacobs 2013 pic a based scale
    Trms = n.ma.masked_where(F['mask'],F['Trms'])*flux_scale
    if FILENAME_LABELS:
        label=file
        ls={}
    else:
        try:
            if F['label'] != 'None':
                print F['label']
                label=str(F['label'])
                ls = {'color':linestyles[label]}
            else:
                label=file
                ls = {}
        except(KeyError):
            label=file
            ls = {}
    print "label = ",label,ls
    if FILENAME_LABELS:
        plot(freqs*1e3,Trms,label=label)
    else:
        plot(freqs*1e3,Trms,label=label,**ls)
    
    if FILENAME_LABELS: legend(loc='upper right')
#    savefig('diff_uv.png')
    grid(which='both')
    ylabel('Brightness Temperature [mK]')
    xlabel('MHz')


#NOISE CALCULATION
hours = 3600 #seconds / hour

kHz = 1e3 #Hz/kHz
Tsky = 450e3 #Sky temp [mK] at 160MHz
f0 = 160.
Trcvr = 100e3 #receiver temp in mK
#nights=20
nights = 92
#nights = 46
channel_width = 492.61 #in kHz
bl_per_column = 1
tint = 43 #integration time in seconds
Npol = 2
Nbls = 7*4.

f = n.linspace(100,200)
for nights in [92]:
    skynoise = Tsky*((f/f0)**(-2.6))/n.sqrt(Npol*nights * tint * (channel_width *kHz) *bl_per_column)
    rcvrnoise = Trcvr/n.sqrt(Npol*nights * tint * (channel_width *kHz) *bl_per_column)
    plot(f,skynoise+rcvrnoise,':k',label='expected noise')
    print "after",nights,"nights sky noise @ 160MHz ", n.interp(160,f,(skynoise+rcvrnoise)/1e3) , "mK"


#noise after averaging bls
#plot(f,(skynoise+rcvrnoise)/Nbls,':k',label='bl avg')


#estimate the fringe rate filtered noise level
t_fr = 1/(15 * w_earth)#|bl in wavelength| * w_earth (see Eq 3 in Psa32_maxred)
intfactor = t_fr/tint
print 'fringe rate filter integration time = [s]',t_fr
plot(f,(skynoise+rcvrnoise)/n.sqrt(intfactor),':k',label='expected noise post fr filter')

#noise level of bl average

#plot(f,(skynoise+rcvrnoise)/n.sqrt(intfactor)/Nbls,':k',label='fr filter + bl avg')


#CLEANING UP THE PLOT
#legend(loc="upper right",prop={'size':12})
ylim([0,TMAX])
xlim([110,190])
if False: #put channels on the top axis
    twiny()
    xlim([0,203])
    xlabel('channel [500kHz]')
elif True: # put redshifts on the top axis
    xlims = ax1.get_xlim()
    ax2 = ax1.twiny()
    x_z = n.array(list(set(n.round(f212z(n.linspace(ax1.get_xlim()[0],ax1.get_xlim()[1])*1e6))))+[7,6])
    x_zfreqs = f21/(1+x_z)/1e6
    ax2.set_xticks(x_zfreqs)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticklabels(x_z)
    ax2.set_xlabel('redshift') 
ylim([TMIN,TMAX])
sca(ax1)
grid(which='both')
if False: #a dumb extra section to plot the noise ratio, change noise model to match the input uv_rms
    noise_model = n.interp(freqs*1e3,f,skynoise+rcvrnoise)/n.sqrt(intfactor)
    figure()
    for file in files:    
        F = n.load(file)
        freqs = F['freqs']
        Trms = n.ma.masked_where(F['mask'],F['Trms'])
        try:
            if F['label'] != 'None':
                label=str(F['label'])
                ls = {'color':linestyles[label]}
            else:
                label=file
                ls = {}
        except(KeyError):
            label=file
            ls = {}
        print "label = ",label,ls
        #if label.endswith('fr'): continue
        plot(freqs*1e3,Trms/noise_model,label=label,**ls)
        
    #    savefig('diff_uv.png')
        grid()
        ylabel('T_measure/T_theory')
        xlabel('MHz')
        ylim([0,4])

show()
