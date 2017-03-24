#!/usr/bin/env python
import matplotlib as mpl
#mpl.rcParams['font.size'] = 12
#mpl.rcParams['xtick.labelsize']= 18
#mpl.rcParams['legend.fontsize'] = 14
#mpl.rcParams['lines.markersize']= 12
#
#mpl.rcParams['legend.numpoints']  = 1
#mpl.rcParams['legend.frameon'] = False
#mpl.rcParams['figure.dpi'] = 500
#mpl.rcParams['savefig.dpi'] = 500
#mpl.rcParams['savefig.format'] ='png'
import numpy as np, aipy, sys, os, optparse
import py21cmsense as py21cm
import capo, capo.eor_results, glob
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from IPython import embed

o = optparse.OptionParser()

o.set_usage('plot_noise_curves.py [options]')
o.set_description(__doc__)
o.add_option('--plot',action='store_true',
    help='output plot before saving')
o.add_option('--noisefiles',type='string',
    help='supply 21cmSense files to plot sensitivity')
o.add_option('--outfile', type='string',
    help='give custom output file name')
o.add_option('--ratio', action='store_true',
    help='Plots ratio of pC to pI under curve')
opts,args = o.parse_args(sys.argv[1:])


data={
    'pI': [x for x in args if 'pI' in x],
    'pC': [x for x in args if 'pI' not in x]
        }

colors = {
        "pI" : "blue",
        "pC" : "green"
        }

z,_,_,_ = capo.eor_results.get_k3pk_from_npz( data['pI'] )
try: Nzs = len(z)
except: Nzs=1
if opts.ratio:
    fig = plt.figure()
    gs = gridspec.GridSpec(3,Nzs)
    gs.update(hspace=0.0, wspace=0.2)
    ax1 = [plt.subplot( gs[:-1,i]) for i in range(Nzs)]
    ax2 = [plt.subplot( gs[-1,i] , sharex=ax1[i]) for i in range(Nzs)]

    fig2 = plt.figure()
    gs = gridspec.GridSpec(3,Nzs)
    gs.update(hspace=0.0, wspace=0.2)
    ax3 = [plt.subplot( gs[:-1,i]) for i in range(Nzs)]
    ax4 = [plt.subplot( gs[-1,i] , sharex=ax1[i]) for i in range(Nzs)]
else:
    fig = plt.figure()
    gs = gridspec.GridSpec(1,Nzs)
    gs.update(hspace=0.0, wspace=0.2)
    ax1 = [plt.subplot( gs[:,i]) for i in range(Nzs)]

    fig2 = plt.figure()
    gs = gridspec.GridSpec(1,Nzs)
    gs.update(hspace=0.0, wspace=0.2)
    ax3 = [plt.subplot( gs[:,i]) for i in range(Nzs)]

noise_freqs, noise_ks, noises = py21cm.load_noise_files(
        glob.glob(opts.noisefiles),polyfit_deg=3
        )
POBER_NOISE = py21cm.noise_interp2d(noise_freqs,noise_ks,noises,interp_kind='linear')

for key in data:
    z,ks,k3pk,k3err = capo.eor_results.get_k3pk_from_npz( data[key] )
    _,kpls,pk,pkerr = capo.eor_results.get_pk_from_npz( data[key] )
    fqs = capo.pspec.z2f(z)*1e3 ##put freqs in Mhz
    for i, redshift in enumerate(z):

        ax1[i].plot(ks[i], np.abs(k3pk[i]) + k3err[i],
                '--', color=colors[key],label=key)

        ax1[i].plot(ks[i], 2*POBER_NOISE(fqs[i],ks[i]), 'k-')

        ax1[i].set_yscale('log')
        if i > 0:
            ax1[i].get_shared_y_axes().join(ax1[0],ax1[i])
        ax1[i].set_title('z = {0:.2f}'.format(redshift))

        if i==0:
            ax1[i].set_ylabel('$\\frac{k^{3}}{2\pi^{2}} P(k) [mK]^{2}$')
        if opts.ratio: plt.setp( ax1[i].get_xticklabels(), visible=False)
        ax1[i].grid(True)

        d2_n = 2*POBER_NOISE(fqs[i],ks[i])
        pk_n = 2*np.pi**2/(ks[i]**3) * d2_n

        ax3[i].plot(kpls[i], np.abs(pk[i]) + pkerr[i],
                '--', color=colors[key],label=key)
        # print len(kpls[i]), len(pk_n)
        if len(kpls[i]) > len(pk_n):
            truncate = len(kpls[i]) - len(pk_n)
            plt_k = kpls[i][truncate:]
            ax3[i].plot(plt_k, pk_n, '-', color='black')
            ax3[i].plot(-plt_k, pk_n, '-', color='black')
        else: ax3[i].plot(kpls[i], pk_n, '-',color='black')
        ax3[i].set_yscale('log')

        ax3[i].set_title('z = {0:.2f}'.format(redshift))
        if i > 0:
            ax3[i].get_shared_y_axes().join(ax3[0],ax3[i])

        if i==0:
            ax3[i].set_ylabel('$ P(k) \\frac{[mK]^{2}}{(hMpc^{-1})^{3}}$')
        if opts.ratio: plt.setp( ax3[i].get_xticklabels(), visible=False)
        ax3[i].grid(True)

if opts.ratio:
    print 'z\tpC/pI\tstd'
    for i,redshift in enumerate(z):
        z_pI,ks_pI,k3pk_pI,k3err_pI = capo.eor_results.get_k3pk_from_npz( data['pI'] )
        z_pC,ks_pC,k3pk_pC,k3err_pC = capo.eor_results.get_k3pk_from_npz( data['pC corrected'] )
        fqs = capo.pspec.z2f(z_pI)*1e3 ##put freqs in Mhz
        ratio = ( np.abs(k3pk_pC[i]) + k3err_pC[i])/(np.abs(k3pk_pI[i]) + k3err_pI[i])
        ax2[i].plot(ks_pI[i], ratio,'-')

        print '{0:.2f}\t{1:.3f}\t{2:.3f}'.format( redshift, np.mean(ratio), np.std(ratio))

        if i==0:
            ax2[i].set_ylabel('pC/pI')
        ax2[i].set_xlabel('$k [hMpc^{-1}]$')
        #ax2[i].set_yscale('log')
        nbins = len(ax2[i].get_yticklabels()) # added
        ax2[i].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
        if i==0:
            nbins = len(ax2[0].get_xticklabels()) # added
            ax2[0].xaxis.set_major_locator(MaxNLocator(nbins=nbins-3))
        ax2[i].grid(True)

        _,kpls_pI,pk_pI,pkerr_pI = capo.eor_results.get_pk_from_npz( data['pI'] )
        _,kpls_pC,pk_pC,pkerr_pC = capo.eor_results.get_pk_from_npz( data['pC corrected'] )
        ratio2 = ( np.abs(pk_pC[i]) + pkerr_pC[i])/(np.abs(pk_pI[i]) + pkerr_pI[i])
        ax4[i].plot(kpls_pI[i], ratio2,'-')

        print '{0:.2f}\t{1:.3f}\t{2:.3f}'.format( redshift, np.mean(ratio), np.std(ratio))

        if i==0:
            ax4[i].set_ylabel('pC/pI')
        ax4[i].set_xlabel('$k [hMpc^{-1}]$')
        #ax2[i].set_yscale('log')
        nbins = len(ax4[i].get_yticklabels()) # added
        ax4[i].yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
        if i==0:
            nbins = len(ax4[0].get_xticklabels()) # added
            ax4[0].xaxis.set_major_locator(MaxNLocator(nbins=nbins-3))
        ax4[i].grid(True)
handles, labels = ax1[-1].get_legend_handles_labels()
ax1[-1].legend(reversed(handles), reversed(labels), loc='lower right')  # reverse
handles, labels = ax3[-1].get_legend_handles_labels()
ax3[-1].legend(reversed(handles), reversed(labels), loc='lower right')  # reverse to keep order consistent
# ax1[-1].legend(loc='lower right')
# ax3[-1].legend(loc='lower right')
#fig.legend(lines,labels,'best')
#gs.tight_layout(fig)
if opts.outfile is not None:
    fig.savefig(opts.outfile+'.png',format='png')
    fig2.savefig(opts.outfile+'_pk.png',format='png')
else:
    fig.savefig('noise_curve_plot.png',format='png')
    fig2.savefig('noise_curve_plot_pk.png',format='png')
if opts.plot: plt.show()
# embed()
plt.close()
