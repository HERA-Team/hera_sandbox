#!/usr/bin/env python
"""Create simple 2 sigma upper limit plot."""

import matplotlib as mpl
import numpy as np
import aipy
import sys
import os
import optparse
import py21cmsense as py21cm
import capo
import capo.eor_results
import glob
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from IPython import embed
# mpl.rcParams['font.size'] = 12
# mpl.rcParams['xtick.labelsize']= 18
# mpl.rcParams['legend.fontsize'] = 14
# mpl.rcParams['lines.markersize']= 12
#
# mpl.rcParams['legend.numpoints']  = 1
# mpl.rcParams['legend.frameon'] = False
# mpl.rcParams['figure.dpi'] = 500
# mpl.rcParams['savefig.dpi'] = 500
# mpl.rcParams['savefig.format'] ='png'

o = optparse.OptionParser()

o.set_usage('plot_noise_curves.py [options]')
o.set_description(__doc__)
o.add_option('--plot', action='store_true',
             help='output plot before saving')
o.add_option('--noisefiles', type='string',
             help='supply 21cmSense files to plot sensitivity')
o.add_option('--outfile', type='string',
             help='give custom output file name')
o.add_option('--ratio', action='store_true',
             help='Plots ratio of pC to pI under curve')
o.add_option('--psa32', action='store_true',
             help='Plots psa32 results to compare')
o.add_option('--psa32_noise', type='string',
             help='Plots noise from psa32')
opts, args = o.parse_args(sys.argv[1:])


data = {
    'unweighted': [x for x in args if 'pI' in x],
    'weighted': [x for x in args if 'pI' not in x]
        }

# colors = {
#         "unweighted": "blue",
#         "weighted": "green"
#         }

z, _, _, _ = capo.eor_results.get_k3pk_from_npz(data['unweighted'])
z_ref, z_ind = np.unique(z, return_index=True)
# reverse because unique orders ascending
z_ref = z_ref[::-1]
z_ind = z_ind[::-1]
# Trim off uncessary unweighted spectra
data['unweighted'] = np.take(data['unweighted'], z_ind, axis=0).tolist()

try:
    Nzs = len(z_ref)
except:
    Nzs = 1
if opts.ratio:
    fig = plt.figure()
    gs = gridspec.GridSpec(3, Nzs)
    gs.update(hspace=0.0, wspace=0.2)
    ax1 = [plt.subplot(gs[:-1, i]) for i in range(Nzs)]
    ax2 = [plt.subplot(gs[-1, i], sharex=ax1[gs_ind]) for i in range(Nzs)]

    fig2 = plt.figure()
    gs = gridspec.GridSpec(3, Nzs)
    gs.update(hspace=0.0, wspace=0.2)
    ax3 = [plt.subplot(gs[:-1, i]) for i in range(Nzs)]
    ax4 = [plt.subplot(gs[-1, i], sharex=ax1[gs_ind]) for i in range(Nzs)]
else:
    fig = plt.figure()
    gs = gridspec.GridSpec(1, Nzs)
    gs.update(hspace=0.0, wspace=0.2)
    ax1 = [plt.subplot(gs[:, i]) for i in range(Nzs)]

    fig2 = plt.figure()
    gs = gridspec.GridSpec(1, Nzs)
    gs.update(hspace=0.0, wspace=0.2)
    ax3 = [plt.subplot(gs[:, i]) for i in range(Nzs)]

noise_freqs, noise_ks, noises = py21cm.load_noise_files(
        glob.glob(opts.noisefiles), polyfit_deg=3
        )
freqs = capo.pspec.z2f(z_ref) * 1e3
noise_ind = np.array(
                     [np.argmin(abs(np.array(noise_freqs) - fq))
                      for fq in freqs])

noise_freqs = np.take(noise_freqs, noise_ind).tolist()
noise_ks = np.take(noise_ks, noise_ind, axis=0).tolist()
noises = np.take(noises, noise_ind, axis=0).tolist()

# POBER_NOISE = py21cm.noise_interp2d(
#               noise_freqs,noise_ks,noises,interp_kind='linear')

if opts.psa32_noise:
    noise_freqs_32, noise_ks_32, noises_32 = py21cm.load_noise_files(
            glob.glob(opts.psa32_noise), polyfit_deg=3
            )
    # POBER_NOISE_32 = py21cm.noise_interp2d(
#              noise_freqs_32,noise_ks_32,noises_32,interp_kind='linear')


for key in data:
    if not data[key]:
        continue
    z, ks, k3pk, k3err = capo.eor_results.get_k3pk_from_npz(data[key])
    _, kpls, pk, pkerr = capo.eor_results.get_pk_from_npz(data[key])
    fqs = capo.pspec.z2f(z)*1e3  # put freqs in Mhz
    k_max = np.max(ks)
    for i, redshift in enumerate(z):

        if key == 'unweighted':
            label = key
        else:
            prob = data[key][i].split('_')[-1].split('.')[0]
            label = key + ' ' + prob

        gs_ind = np.where(z_ref == redshift)[0].item()
        # special index for gridspec

        ax1[gs_ind].plot(ks[i], k3pk[i] + k3err[i],
                         '--', label=label)

        ax1[gs_ind].set_yscale('log')
        if i > 0:
            ax1[gs_ind].get_shared_y_axes().join(ax1[0], ax1[gs_ind])
        ax1[gs_ind].set_title('z = {0:.2f}'.format(redshift))
        ax1[gs_ind].set_xlabel('$k$ [$h$ Mpc$^{-1}$]')

        if i == 0:
            ax1[gs_ind].set_ylabel('$\\frac{k^{3}}{2\pi^{2}} P(k) [mK]^{2}$')
            ax1[0].set_ylim([1e0, 1e9])
        if opts.ratio:
            plt.setp(ax1[gs_ind].get_xticklabels(), visible=False)
        ax1[gs_ind].grid(True)
        ax1[gs_ind].set_xlim(0, k_max * 1.01)

        if Nzs == 1:
            if i == 0:
                ax1[gs_ind].plot(noise_ks[0], noises[0], 'k-')
            else:
                ax1[gs_ind].plot(noise_ks[0], noises[0], 'k-')  # * 3887./2022.
            d2_n = noises[0]
            pk_n = 2*np.pi**2/(np.array(noise_ks[0])**3)*d2_n
            if len(kpls[0]) > len(pk_n):
                ax3[gs_ind].plot(noise_ks[0], pk_n, '-', color='black')
                ax3[gs_ind].plot(-noise_ks[0], pk_n, '-', color='black')
            else:
                ax3[gs_ind].plot(noise_ks[0], pk_n, '-', color='black')
        else:
            ax1[gs_ind].plot(noise_ks[gs_ind], noises[gs_ind], 'k-')
            d2_n = noises[gs_ind]
            pk_n = 2*np.pi**2/(np.array(noise_ks[gs_ind])**3)*d2_n
            if len(kpls[gs_ind]) > len(pk_n):
                ax3[gs_ind].plot(noise_ks[gs_ind], pk_n, '-', color='black')
                ax3[gs_ind].plot(-noise_ks[gs_ind], pk_n, '-', color='black')
            else:
                ax3[gs_ind].plot(noise_ks[gs_ind], pk_n, '-', color='black')

        ax3[gs_ind].plot(kpls[i], pk[i] + pkerr[i],
                         '--', label=label)
        # print len(kpls[i]), len(pk_n)

        # add optional psa32 data
        if key is data.keys()[-1]:
            if opts.psa32:
                PSA32_z, PSA32_pspec = capo.eor_results.z_slice(
                            redshift, capo.eor_results.PAPER_32_all())
                # print "loading GMRT 2014 data near redshift:",GMRT_z
                if np.abs(PSA32_z - redshift) < 0.5:
                    ax1[gs_ind].plot(PSA32_pspec[:, 0], PSA32_pspec[:, 2],
                                     '--', label='psa32')

                    psa32_k = PSA32_pspec[:, 0]
                    psa32_pk = PSA32_pspec[:, 2]*2*np.pi**2/psa32_k**3
                    ax3[gs_ind].plot(PSA32_pspec[:, 0], psa32_pk,
                                     '--', label='psa32')
            if opts.psa32_noise:
                ax1[gs_ind].plot(noise_ks_32[i], noises_32[i], 'k-.')

                d2_n_32 = noises_32[i]
                pk_n_32 = 2*np.pi**2/(np.array(noise_ks_32[i])**3) * d2_n_32

                if len(kpls[i]) > len(pk_n_32):
                    ax3[gs_ind].plot(noise_ks_32[i],
                                     pk_n_32, '-.', color='black')
                    ax3[gs_ind].plot(-noise_ks_32[i],
                                     pk_n_32, '-.', color='black')
                else:
                    ax3[gs_ind].plot(noise_ks_32[i],
                                     pk_n_32, '-.', color='black')

        ax3[gs_ind].set_yscale('log')

        ax3[gs_ind].set_title('z = {0:.2f}'.format(redshift))
        if i > 0:
            ax3[gs_ind].get_shared_y_axes().join(ax3[0], ax3[gs_ind])

        if i == 0:
            ax3[gs_ind].set_ylabel('$ P(k) \\frac{[mK]^{2}}{(hMpc^{-1})^{3}}$')
        if opts.ratio:
            plt.setp(ax3[gs_ind].get_xticklabels(), visible=False)
        ax3[gs_ind].grid(True)
        ax3[gs_ind].set_xlim(0, k_max * 1.01)


if opts.ratio:
    print 'z\tpC/pI\tstd'
    for i, redshift in enumerate(z):
        z_pI, ks_pI, k3pk_pI, k3err_pI = capo.eor_results.get_k3pk_from_npz(
                                        data['pI'])
        z_pC, ks_pC, k3pk_pC, k3err_pC = capo.eor_results.get_k3pk_from_npz(
                                        data['pC corrected'])
        fqs = capo.pspec.z2f(z_pI)*1e3  # put freqs in Mhz
        ratio = (k3pk_pC[i] + k3err_pC[i])/(
                    k3pk_pI[i] + k3err_pI[i])
        ax2[gs_ind].plot(ks_pI[i], ratio, '-')

        print '{0:.2f}\t{1:.3f}\t{2:.3f}'.format(
                        redshift, np.mean(ratio), np.std(ratio))

        if i == 0:
            ax2[gs_ind].set_ylabel('pC/pI')
        ax2[gs_ind].set_xlabel('$k [hMpc^{-1}]$')
        # ax2[gs_ind].set_yscale('log')
        nbins = len(ax2[gs_ind].get_yticklabels())  # added
        ax2[gs_ind].yaxis.set_major_locator(
                                    MaxNLocator(nbins=nbins, prune='upper'))
        if i == 0:
            nbins = len(ax2[0].get_xticklabels())  # added
            ax2[0].xaxis.set_major_locator(MaxNLocator(nbins=nbins-3))
        ax2[gs_ind].grid(True)
        ax2[gs_ind].set_xlim(0, k_max * 1.01)

        _, kpls_pI, pk_pI, pkerr_pI = capo.eor_results.get_pk_from_npz(
                                    data['pI'])
        _, kpls_pC, pk_pC, pkerr_pC = capo.eor_results.get_pk_from_npz(
                                    data['pC corrected'])
        ratio2 = (pk_pC[i] + pkerr_pC[i])/(
                        pk_pI[i] + pkerr_pI[i])
        ax4[gs_ind].plot(kpls_pI[i], ratio2, '-')

        print '{0:.2f}\t{1:.3f}\t{2:.3f}'.format(
                    redshift, np.mean(ratio), np.std(ratio))

        if i == 0:
            ax4[gs_ind].set_ylabel('pC/pI')
        ax4[gs_ind].set_xlabel('$k [hMpc^{-1}]$')
        # ax2[i].set_yscale('log')
        nbins = len(ax4[gs_ind].get_yticklabels())  # added
        ax4[gs_ind].yaxis.set_major_locator(
                                MaxNLocator(nbins=nbins, prune='upper'))
        if i == 0:
            nbins = len(ax4[0].get_xticklabels())  # added
            ax4[0].xaxis.set_major_locator(MaxNLocator(nbins=nbins-3))
        ax4[gs_ind].grid(True)
        ax4[gs_ind].set_xlim(0, k_max * 1.01)

handles, labels = ax1[-1].get_legend_handles_labels()
ax1[-1].legend(reversed(handles), reversed(labels),
               loc='lower right', numpoints=1)  # reverse
handles, labels = ax3[-1].get_legend_handles_labels()
ax3[-1].legend(reversed(handles), reversed(labels),
               loc='upper right', numpoints=1)
# reverse to keep order consistent
# ax1[-1].legend(loc='lower right')
# ax3[-1].legend(loc='lower right')
# fig.legend(lines,labels,'best')
# gs.tight_layout(fig)
if opts.outfile is not None:
    fig.savefig(opts.outfile+'.png', format='png')
    fig2.savefig(opts.outfile+'_pk.png', format='png')
else:
    fig.savefig('noise_curve_plot.png', format='png')
    fig2.savefig('noise_curve_plot_pk.png', format='png')
if opts.plot:
    plt.show(block=1)
# embed()
# plt.close()
