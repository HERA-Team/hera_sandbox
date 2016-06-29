#! /usr/bin/env python
"""
plot the output frf_diagnose.py pickles
"""
import matplotlib as mpl
import numpy as np

#mpl.rcParams['font.size']  = 12
mpl.rcParams['legend.numpoints']  = 1
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.fontsize'] = 12
#mpl.rcParams['figure.dpi'] = 300
#mpl.rcParams['savefig.dpi'] = 300
#mpl.rcParams['savefig.format'] ='png'
#mpl.rcParams['lines.markeredgewidth'] = 0
mpl.rcParams['lines.linewidth'] = 4
mpl.rcParams['figure.figsize'] = (12,5)
import pickle
import matplotlib.cm as cmx
import matplotlib.colors as colors
import numpy as np, matplotlib.pyplot as p, glob, optparse, sys

o = optparse.OptionParser()
o.set_usage('plot_frf.py [options]')
o.set_description(__doc__)
o.add_option('--plot',action='store_true',
    help='outputs plots before saving')
o.add_option('--chan',type=int,default=101,
    help='channel index to plot')
opts,args = o.parse_args(sys.argv[1:])

def get_colors(N):
    '''Returns function with N unique colors'''
    norm=colors.Normalize(vmin=0,vmax=N-1)
    scal_map=cmx.ScalarMappable(norm=norm,cmap='brg')
    def map_index_to_rgb(index):
        return scal_map.to_rgba(index)
    return map_index_to_rgb

#mychan = opts.chan
mysep = '0,1'
fig, ax = p.subplots(1)
fig2,ax2 = p.subplots(1)
#ax.set_title('Fitted Fringe Rate Profile')
ax.set_xlabel('Fringe Rate [mili Hz]')
ax.set_xlim([-.7,1.5])
ax.set_ylim([0,1])

ax2.set_xlabel('time [s]')
ax2.set_xlim([-10000,10000])
#ax2.set_ylim([-1,1])
#ax2.set_yscale('log')
ax2.set_xlabel('time [s]')

num_files = len(args) + 2

cmap = get_colors(num_files)

for cnt,filename in enumerate(args):
    print("plotting:{filename}".format(filename=filename))
    F = pickle.load(open(filename,'r'))

    mychan = F['chan']
    frp_freqs = F['freqs']
    frp = F['frp']
    frp_fit = F['frps'][mychan]
    frp_fit /= np.max(frp_fit)
    fir = F['firs']
    timebins = F['timebins']
    fir_single = fir[mychan]
    #factor = np.sqrt(np.sum(np.abs(fir_single)**2))
    #fir_single /= np.abs(fir_single).max()
    #fir_single /= factor


    fir_fft = np.fft.ifft(np.fft.ifftshift(frp_fit))
    fir_fft = np.fft.fftshift(fir_fft)
    fir_fft /=  np.abs(fir_fft).max()

    lines = ax.plot( frp_freqs * 1e3, frp, label =filename,color=cmap(cnt),alpha=.5)
    color = lines[0].get_color()
    ax.plot(frp_freqs*1e3,frp_fit,'--',color=cmap(cnt),alpha=.5)
    ax2.plot(timebins,fir_single.real,'-',color=cmap(cnt),label=filename,alpha=.5)
    #ax2.plot(timebins,fir_single.imag,':',color=color)
    #ax2.plot(timebins,fir_fft,'--',color=color)
    #ax2.plot(timebins,fir_fft.imag,'-.',color=color)
ax2.legend(loc='upper right')
ax.legend(loc='best')
fig.savefig('frf_comparison.png',format='png')
fig2.savefig('fir_comparison_complex.png')


if opts.plot:
    p.show()

