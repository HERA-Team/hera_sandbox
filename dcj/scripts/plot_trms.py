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
TMAX=40 #maximum plot in mK

w_earth = 2*pi*1.16e-5 #sidereal rotation rate of earth in rad/s
matplotlib.rcParams.update({'font.size':18})


#LOAD AND PLOT THE output from uv_rms.py
files = sys.argv[1:]
linestyles = {'bl':'b','bl,fr':'b','df':'g','df,fr':'g'}
ax1 = subplot(111)
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
    plot(freqs*1e3,Trms,label=label,**ls)

#    savefig('diff_uv.png')
    grid()
    ylabel('mK')
    xlabel('MHz')


#NOISE CALCULATION
hours = 3600 #seconds / hour

kHz = 1e3 #Hz/kHz
Tsky = 450e3 #Sky temp [mK] at 150MHz
Trcvr = 100e3 #receiver temp in mK
nights = 92
channel_width = 492.61 #in kHz
bl_per_column = 1
tint = 39 #integration time in seconds
Npol = 2

f = n.linspace(100,200)
skynoise = Tsky/n.sqrt(Npol*nights * tint * (channel_width *kHz) *bl_per_column)*(f/150.)**(-2.6)
rcvrnoise = Trcvr/n.sqrt(Npol*nights * tint * (channel_width *kHz) *bl_per_column)
plot(f,skynoise+rcvrnoise,':k',label='expected noise')

#estimate the fringe rate filtered noise level
t_fr = 1/(15 * w_earth)#|bl in wavelength| * w_earth (see Eq 3 in Psa32_maxred)
intfactor = t_fr/tint
print 'fringe rate filter integration time = [s]',t_fr
plot(f,(skynoise+rcvrnoise)/n.sqrt(intfactor),':k',label='expected noise post fr filter')

#CLEANING UP THE PLOT
#legend(loc="upper right",prop={'size':12})
ylim([0,TMAX])
xlim([110,190])
grid()
if False: #put channels on the top axis
    twiny()
    xlim([0,203])
    xlabel('channel [500kHz]')
else: # put redshifts on the top axis
    xlims = ax1.get_xlim()
    ax2 = ax1.twiny()
    x_z = n.array(list(set(n.round(f212z(n.linspace(ax1.get_xlim()[0],ax1.get_xlim()[1])*1e6))))+[7,6])
    x_zfreqs = f21/(1+x_z)/1e6
    ax2.set_xticks(x_zfreqs)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticklabels(x_z)
    ax2.set_xlabel('redshift') 

if True: #a dumb extra section to plot the noise ratio
    noise_model = n.interp(freqs*1e3,f,skynoise+rcvrnoise)
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
        if label.endswith('fr'): continue
        plot(freqs*1e3,Trms/noise_model,label=label,**ls)
        
    #    savefig('diff_uv.png')
        grid()
        ylabel('T_measure/T_theory')
        xlabel('MHz')

show()
