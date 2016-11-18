#! /usr/bin/env python

import aipy, numpy as np, pylab as plt, capo
import sys

info,data,flgs = capo.miriad.read_files(sys.argv[1:], antstr='auto', polstr='xx,yy')
POL = data.values()[0].keys()[0]
CHUNK = 256
BINS = 24
colors = 'bgrcmy'
N = np.ceil(np.sqrt(len(data)))
M = np.ceil(len(data) / float(N))
bins = np.logspace(-2,4,BINS,base=2.)
ants = [i for (i,_) in data]; ants.sort()
for cnt,i in enumerate(ants):
    ax = plt.subplot(N,M,cnt+1)
    for j,ch in enumerate(xrange(0,1024,CHUNK)):
        d = data[(i,i)][POL][:,ch:ch+CHUNK].flatten()
        h,b = np.histogram(np.sqrt(d/2),bins)
        h = 10**np.log10(h+.1)
        b = 0.5*(b[1:] + b[:-1])
        ax.fill_between(np.log2(b), h, .1, where=h>.1, color=colors[j], alpha=.5)
    bounds = np.where(bins < 2**0, d.size, np.where(bins > 2**2, d.size, 0))
    ax.fill_between(np.log2(bins), bounds, .1, where=bounds>.1, color='black', alpha=.5)
    ax.set_yscale('log')
    plt.xlim(-2,3)
    plt.ylim(d.size/1e2, d.size)
    plt.title(str(i)+POL)
    ax.get_yaxis().set_visible(False)
    if cnt < (N-1)*M:
        ax.get_xaxis().set_ticklabels([])
    else: plt.xlabel(r'$V_{\rm rms}$ [bits]')
    plt.grid()
plt.subplots_adjust(wspace=.05, hspace=.4)
plt.show()
    
