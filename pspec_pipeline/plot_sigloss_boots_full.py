#! /usr/bin/env python
"""Takes injections and creates output distributions."""
import matplotlib
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z, dk_du
import py21cmsense as py21cm
import sigloss_functions as sf
import glob
import ipdb
import sys
from IPython import embed


def errorbars(data, axis=1, per=95):
    """Return errorbars as percnitles of distributions."""
    mean = n.percentile(data, 50, axis=axis)
    lower = mean - n.percentile(data, 50-per/2., axis=axis)
    upper = n.percentile(data, 50+per/2., axis=axis) - mean
    lower, upper = n.std(data, axis=axis), n.std(data, axis=axis)
    return lower, uppers


def get_2dhist(p1, p2, log=True, bins=(300, 300)):
    """Creates2-D histogram with convenient parameters."""
    if log:
        p1, p2 = n.log10(p1), n.log10(p2)
    H, xbins, ybins = n.histogram2d(p1.flatten(), p2.flatten(),
                                    bins=bins)
    xb, yb = (xbins[1:] + xbins[:-1])/2., (ybins[1:] + ybins[:-1])/2.
    if log:
        xb, yb, = 10**(xb), 10**(yb)
    H = H.T
    return xb, yb, H


def gauss(x, u, s):
    """Quick gaussian."""
    return n.exp(-(x-u)**2/(2.*s**2))


def interp_gauss(x, func):
    """Interp range of gaussians."""
    return interp1d(x, func, kind='linear', bounds_error=False, fill_value=0)


NBOOT = 100
t_eff = 69
# GETTING SIGLOSS DATA ###
pCs, pIs, pCvs = [], [], []
pIvs = []
pCs_err, pIs_err = [], []
freq = []
kpls = []
files = glob.glob('inject_*')
nums = [float(f.split('_')[-1]) for f in files]
inds = n.argsort(nums)
nums = [nums[ind] for ind in inds]
files = [files[ind] for ind in inds]
for inject in files:
    print 'Reading', inject
    pspecs = glob.glob(inject + '/pspec_boot*.npz')
    inject = float(inject.split('_')[-1])
    pC_avg, pI_avg, pCv_avg = [], [], []
    pIv_avg = []
    pC_spec, pI_spec = [], []
    _kpls = []
    for pspec in pspecs:
        npz = n.load(pspec)
        try:
            freq = npz['freq']
        except:
            pass
        # Arrange in  (#chan, #times)
        pC, pI, pCv = npz['pk_vs_t'], npz['nocov_vs_t'], npz['pCv']
        pIv = npz['pIv']
        _kpls.append(npz['kpl'])
        # sigloss_sim.py makes pC =pC_(data+eor) - pC_data ###
        # This next line sets pC = p_out ###
        pC_avg.append(pC.real[:, ::t_eff])
        # + pCv.real) #avg over freq and time
        pI_avg.append(pI.real[:, ::t_eff])
        pIv_avg.append(pIv.real[:, ::t_eff])
        pCv_avg.append(pCv.real[:, ::t_eff])  # spectrum
    kpls.append(_kpls)
    pCs.append(pC_avg)
    # should bootstrap all these pC_avg's (20 of them) to get error bars...
    # look at their distribution
    pIs.append(pI_avg)
    pIvs.append(pIv_avg)
    pCvs.append(pCv_avg)  # spectrum

pIs, pCs, pCvs, pIvs = n.array(pIs), n.array(pCs), n.array(pCvs), n.array(pIvs)
kpls = n.array(kpls[0][0][:])
pCvs_pk, pIvs_pk = n.copy(pCvs), n.copy(pIvs)

# Make format: inj, k, bootstrap*time
pIs = n.reshape(pIs, (pIs.shape[0], pIs.shape[2], pIs.shape[1]*pIs.shape[3]))
pCs = n.reshape(pCs, (pCs.shape[0], pCs.shape[2], pCs.shape[1]*pCs.shape[3]))
pCvs = n.reshape(pCvs,
                 (pCvs.shape[0], pCvs.shape[2], pCvs.shape[1]*pCvs.shape[3]))
pIvs = n.reshape(pIvs,
                 (pIvs.shape[0], pIvs.shape[2], pIvs.shape[1]*pIvs.shape[3]))
# P(k) values have k, inj, boots, times
pCvs_pk = n.reshape(pCvs_pk,
                    (pCvs_pk.shape[2], pCvs_pk.shape[0],
                     pCvs_pk.shape[1], pCvs_pk.shape[3]))
pIvs_pk = n.reshape(pIvs_pk,
                    (pIvs_pk.shape[2], pIvs_pk.shape[0],
                     pIvs_pk.shape[1], pIvs_pk.shape[3]))


pIs, pCs, pCvs = (n.ma.masked_invalid(n.array(pIs)),
                  n.ma.masked_invalid(n.array(pCs)),
                  n.ma.masked_invalid(n.array(pCvs)))
pIs_mask, pCs_mask, pCvs_mask = pIs.mask, pCs.mask, pCvs.mask

# Set All negative values to zero
# While not technically correct this is for the contours
pIs, pCs, pCvs = (n.ma.masked_less(pIs, 0),
                  n.ma.masked_less(pCs, 0),
                  n.ma.masked_less(pCvs, 0))

pIs.mask = n.logical_or(pIs.mask, pIs_mask)
pCs.mask = n.logical_or(pCs.mask, pCs_mask)
pCvs.mask = n.logical_or(pCvs.mask, pCvs_mask)
pIs.fill_value = 1e-5
pCs.fill_value = 1e-5
pCvs.fill_value = 1e-5

pIvs = n.ma.masked_invalid(n.array(pIvs))
pIvs_mask = pCvs.mask

pIvs = n.ma.masked_less(pIvs, 0)
pIvs.mask = n.logical_or(pIvs.mask, pIvs_mask)
pIvs.fill_value = 1e-5

pIvs_pk = n.ma.masked_invalid(pIvs_pk)
pCvs_pk = n.ma.masked_invalid(pCvs_pk)
pIvs_pk_mask, pCvs_pk_mask = pIvs_pk.mask, pCvs_pk.mask

pIvs_pk = n.ma.masked_less(pIvs_pk, 0)
pCvs_pk = n.ma.masked_less(pCvs_pk, 0)

pIvs_pk.mask = n.logical_or(pIvs_pk.mask, pIvs_pk_mask)
pCvs_pk.mask = n.logical_or(pCvs_pk.mask, pCvs_pk_mask)
pIvs_pk.fill_value = 1e-5
pCvs_pk.fill_value = 1e-5


try:
    z_bin = f2z(freq)
except:
    print 'frequency not found in boots. Searching in pspec.npz'
    f_name = ('/home/mkolopanis/psa64/sigloss_verification/'
              'Jul6_noise_3Jy_inttime_44/95_115/I/'
              'pspec_Jul6_noise_3Jy_inttime_44_95_115_I.npz')
    npz = n.load(f_name)  # matt's data
    freq = npz['freq']
    z_bin = f2z(freq)

# load 21cmSense Noise models used to compute Beta
# Beta = (P_noise + P_inj)/P_out

noise_files = ('/home/mkolopanis/psa64/'
               '21cmsense_noise/dish_size_1/*drift_mod*.npz')
n_fqs, n_ks, noise = py21cm.load_noise_files(glob.glob(noise_files))

noise_interp = py21cm.noise_interp2d(n_fqs, n_ks, noise)

kpls_pos = n.concatenate(n.array_split(kpls, [10, 11])[1:])
umag = 30/(3e8/(freq*1e9))
kpr = dk_du(z_bin) * umag
n_k = n.array(n.sqrt(kpls_pos**2 + kpr**2))
d2_n = noise_interp(freq*1e3, n_k)
pk_n = 2*n.pi**2/(n_k**3) * d2_n  # * 3887./2022
p_n = n.median(pk_n)
# p_n = n.max( pk_n )

print 'This script does not Bootstrap'

# take the median over k's
# Now has shape (inj, boots*times)
pIs = n.ma.median(pIs, 1)
pCs = n.ma.median(pCs, 1)
pIvs = n.ma.median(pIvs, 1)
pCvs = n.ma.median(pCvs, 1)
pCvs_pk = n.ma.median(pCvs_pk, 1)
pIvs_pk = n.ma.median(pIvs_pk, 1)

# all arrays should be of the form (num ks, num injects) now median over proper
# dimension
pIs_full = pIs.copy()
pCs_full = pCs.copy()
pCvs_full = pCvs.copy()
pIvs_full = pIvs.copy()

# Take Mean and std over boots, times for pCvs_pk, pIvs_pk
pCvs_pk = n.ma.mean(pCvs_pk, -1)
pCvs_pk_err = n.ma.std(pCvs_pk, -1)
pCvs_pk = n.ma.mean(pCvs_pk, -1)

pIvs_pk = n.ma.mean(pIvs_pk, -1)
pIvs_pk_err = n.ma.std(pIvs_pk, -1)
pIvs_pk = n.ma.mean(pIvs_pk, -1)

fig = p.figure(figsize=(7, 7))
gs = gridspec.GridSpec(2, 3, width_ratios=[.4, 1.2, .4],
                       height_ratios=[.4, 1])  # used to be 2,2
fig.subplots_adjust(left=.15, top=.95, bottom=.15,
                    wspace=.35, hspace=.15, right=0.95)

# Plot 2
# p.figure()
pklo, pkhi = 1e-4, 1e18
# pklo,pkhi = 1e-10,1e5 #for noise only
ax2 = p.subplot(gs[4])  # used to be 2
# p.loglog(pIs, pCs, 'k.')
# uncomment if no left-hand P(k) plot
p.setp(ax2.get_yticklabels(), visible=False)
ax2.set_yscale('log', nonposy='clip')

# Get 2-D hist for P_in (pIs) and P_out (pCs)

pIs_bins, pCs_bins, hist_pCs = get_2dhist(pIs.filled(), pCs.filled(),
                                          bins=(300, 300), log=True)
pIs_bins, pIvs_bin, hist_pIvs = get_2dhist(pIs.filled(), pIvs.filled(),
                                           bins=(300, 300), log=True)

p.pcolor(pIs_bins, pCs_bins, hist_pCs, alpha=.9)
p.axhline(p_n, color='green', linestyle='--')

p.xscale('log')
ax2.plot([pklo, pkhi], [pklo, pkhi], 'k-')
p.xlim(pklo, pkhi)
p.ylim(pklo, pkhi)
p.xlabel(r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize=14)
p.ylabel(r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize=14)
p.grid()
pkup = n.max(n.abs(pCvs))
pkdn = n.min(n.abs(pCvs))
# p.fill_between([pklo,pkhi],[pkdn,pkdn],[pkup,pkup],
#                facecolor='gray', edgecolor='gray',alpha=.4)
# p.show(block=False)


# Plot 3
ax3 = p.subplot(gs[5])  # used to be 3
p.setp(ax3.get_yticklabels(), visible=False)
# p.loglog(n.clip(pIs/pCs - 1, 1e-3,n.Inf), pCs, 'k.')
# p.loglog(n.abs(pIs/pCs - 1), pCs, 'k.')

pIs_pn_pCs_bins, pCs_bins, ratio_hist = get_2dhist(
                                    (pIs.filled() + p_n)/pCs.filled(),
                                    pCs.filled(), bins=(300, 300))
p.pcolor(pIs_pn_pCs_bins, pCs_bins, ratio_hist, alpha=.9)

# print "pI/pC : ", pIs/pCs - 1
# print "pC: ", pCs
ax3.set_xscale('log')
ax3.set_yscale('log', nonposy='clip')
p.ylim(pklo, pkhi)
p.xlim(1e-1, 1e4)
# p.xlim(1,10)
p.grid()
p.xlabel(r'$P_{\rm in}/P_{\rm out}$', fontsize=14)
p.fill_between([1e-3, 1e8], [pkdn, pkdn], [pkup, pkup],
               facecolor='gray', edgecolor='gray', alpha=.4)
sig_factors = []
# sig_factors.append(sig_factor_interp(pkup))

# Plot 1
ax0 = p.subplot(gs[1])  # used to be 0
p.setp(ax0.get_xticklabels(), visible=False)

pIs_bins, pCs_pn_pIs_bins, ratio2_hist = get_2dhist(
                            pIs.filled(), pCs.filled()/(pIs.filled() + p_n),
                            bins=(300, 300))

p.pcolor(pIs_bins, pCs_pn_pIs_bins, ratio2_hist, alpha=.9)
ax0.set_xscale('log')
p.xlim(pklo, pkhi)
# p.ylim(3e-2,3)
p.ylim(1e-5, 1e0)
p.yscale('log', nonposy='clip')
p.grid()
p.ylabel(r'$P_{\rm out}/P_{\rm in}$', fontsize=14)

p.title('z = {0:.2f}'.format(z_bin))

# P(k) plot on left side
ax4 = p.subplot(gs[3])
p.setp(ax4.get_yticklabels(), visible=True)
# ax4.set_xscale('log')
ax4.set_yscale('log', nonposy='clip')
p.ylim(pklo, pkhi)

# kpls_bin, pCvs_pk_bin, hist_pCvs = get_2dhist(
#                                 n.tile(kpls, pCvs_pk.shape[1]),
#                                 n.log10(pCvs_pk.filled()),
#                                 log=False)
# kpls_bin, pIvs_pk_bin, hist_pIvs = get_2dhist(
#                                 n.tile(kpls, pIvs_pk.shape[1]),
#                                 n.log10(pIvs_pk.filled()),
#                                 log=False)
# pIvs_pk_bin = 10**(pIvs_pk_bin)
# pCvs_pk_bin = 10**(pCvs_pk_bin)
#
# p.pcolor(kpls_bin, pCvs_pk_bin, hist_pCvs)
# p.contour(kpls_bin, pIvs_pk_bin, hist_pIvs)
p.errorbar(kpls, pCvs_pk, pCvs_pk_err, fmt='b.')
p.errorbar(kpls, pIvs_pk, pIvs_pk_err, fmt='.', color='gray', alpha=.4)
ax0.set_xscale('log')
# p.errorbar(kpls, n.abs(pCvs), yerr=2*pCvs_err, fmt='k.', capsize=0)
p.grid()

p.savefig('sigloss_distribution.png', format='png')

p.show()
