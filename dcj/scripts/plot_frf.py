#! /usr/bin/env python
"""
plot the output frf_diagnose.py pickles
"""
import matplotlib as mpl
import numpy as np

mpl.rcParams['font.size']  = 12
mpl.rcParams['legend.numpoints']  = 1
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['savefig.format'] ='png'
mpl.rcParams['lines.markeredgewidth'] = 0
mpl.rcParams['lines.markersize'] = 7
import pickle
import numpy as np, matplotlib.pyplot as p, glob, optparse, sys

o = optparse.OptionParser()
o.set_usage('plot_frf.py [options]')
o.set_description(__doc__)
o.add_option('--plot',action='store_true',
    help='outputs plots before saving')
o.add_option('--chan',type=int,default=101,
    help='channel index to plot')
opts,args = o.parse_args(sys.argv[1:])

mychan = opts.chan
mysep = '0,1'
fig, ax = p.subplots(1)
fig2,ax2 = p.subplots(1)
#ax.set_title('Fitted Fringe Rate Profile')
ax.set_xlabel('Fringe Rate [mili Hz]')
ax.set_xlim([-.7,1.5])
ax.set_ylim([0,1])

ax2.set_xlabel('time [ns]')
ax2.set_xlim([-1000,1500])
ax2.set_ylim([0.001,1])
ax2.set_yscale('log')
ax2.set_xlabel('time [ns]')

for filename in args:
    print("plotting:{filename}".format(filename=filename))
    F = pickle.load(open(filename,'r'))
    
    frp_freqs = F['freqs']
    frp = F['frp']
    frp_fit = F[mysep][mychan]
    frp_fit /= np.max(frp_fit)
    fir = F['firs']
    timebins = F['timebins']
    fir_single = F[mysep][mychan]
    fir_single /= np.abs(fir_single).max()

    
    lines = ax.plot( frp_freqs * 1e3, frp, label =filename)
    color = lines[0].get_color()
    ax.plot(frp_freqs*1e3,frp_fit,'--',color=color)
    ax2.plot(timebins,fir_single.real,color=color,label=filename)
    ax2.plot(timebins,fir_single.imag,':',color=color)
ax2.legend(loc='best')  
ax.legend(loc='best')
fig.savefig('frf_comparison.png',format='png')
fig2.savefig('fir_comparison_complex.png')


if opts.plot:
    p.show()

