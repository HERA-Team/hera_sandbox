#! /usr/bin/env python
"""Takes Injects and outputs errorbars."""
import matplotlib
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z, dk_du
import py21cmsense as py21cm
#import sigloss_functions as sf
import glob
import ipdb
import sys,numpy as np
from astropy.io import ascii


def errorbars(data, axis=1, per=95):
    """Find upper lower percentiles as errorbars."""
    mean = n.percentile(data, 50, axis=axis)
    lower = mean - n.percentile(data, 50-per/2., axis=axis)
    upper = n.percentile(data, 50+per/2., axis=axis) - mean
    lower, upper = n.std(data, axis=axis), n.std(data, axis=axis)
    return lower, upper


def gauss(x, u, s):
    """"Quick Gaussian."""
    return n.exp(-(x-u)**2/(2.*s**2))


def interp_gauss(x, func):
    """Interpolate Gaussians."""
    return interp1d(x, func, kind='linear', bounds_error=False, fill_value=0)


NBOOT = 100

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
        pC_avg.append(pC.real)  # + pCv.real) #avg over freq and time
        pI_avg.append(pI.real)
        pIv_avg.append(pIv.real)
        pCv_avg.append(pCv.real)  # spectrum
        # pC_spec.append(n.average(pC.real, axis=1))
        # pI_spec.append(n.average(pI.real, axis=1))
        # pC_spec.append(pC.real)
        # pI_spec.append(pI.real)
    kpls.append(_kpls)
    pCs.append(pC_avg)
    # should bootstrap all these pC_avg's (20 of them) to get error bars...
    # look at their distribution
    pIs.append(pI_avg)
    pIvs.append(pIv_avg)
    pCvs.append(pCv_avg)  # spectrum
    # pC_spec = n.average(pC_spec, axis=0)
    # print inject, pspec, pC_avg, pI_avg
    # p.figure(1)
    # p.subplot(211); p.loglog([pI_avg], [pC_avg], 'k.')
    # p.subplot(212); p.loglog([pC_avg], [pI_avg/pC_avg], 'k.')

pIs, pCs, pCvs = (n.ma.masked_invalid(n.array(pIs)),
                  n.ma.masked_invalid(n.array(pCs)),
                  n.ma.masked_invalid(n.array(pCvs)))
pIvs = n.ma.masked_invalid(n.array(pIvs))
pCvs_pk, pIvs_pk = n.copy(pCvs), n.copy(pIvs)
kpls = n.array(kpls[0][0][:])
# inj,k,bootstrap,time
pIs = n.reshape(pIs, (pIs.shape[0], pIs.shape[2], pIs.shape[1], pIs.shape[3]))
pCs = n.reshape(pCs, (pCs.shape[0], pCs.shape[2], pCs.shape[1], pCs.shape[3]))
pCvs = n.reshape(pCvs,
                 (pCvs.shape[0], pCvs.shape[2], pCvs.shape[1], pCvs.shape[3]))
pIvs = n.reshape(pIvs,
                 (pIvs.shape[0], pIvs.shape[2], pIvs.shape[1], pIvs.shape[3]))
# P(k) values have k, inj, boots, times
pCvs_pk = n.reshape(pCvs_pk,
                    (pCvs_pk.shape[2], pCvs_pk.shape[0], pCvs_pk.shape[1],
                     pCvs_pk.shape[3]))
pIvs_pk = n.reshape(pIvs_pk,
                    (pIvs_pk.shape[2], pIvs_pk.shape[0], pIvs_pk.shape[1],
                     pIvs_pk.shape[3]))

# take the median over time
# pIs = n.median(pIs, -1)
# pCs = n.median(pCs, -1)
# pCvs = n.median(pCvs, -1)
# pIvs_pk = n.median(pIvs_pk, -1)
# pCvs_pk = n.median(pCvs_pk, -1)
# pIvs = n.median(pIvs, -1)

try:
    z_bin = f2z(freq)
except:
    print 'frequency not found in boots. Searching in pspec.npz'
    npz = n.load('/home/mkolopanis/psa64/sigloss_verification/'
                 'Jul6_noise_3Jy_inttime_44/95_115/I/'
                 'pspec_Jul6_noise_3Jy_inttime_44_95_115_I.npz')  # matt's data
    freq = npz['freq']
    z_bin = f2z(freq)

# load 21cmSense Noise models used to compute Beta
# Beta = (P_noise + P_inj)/P_out

n_fqs, n_ks, noise = py21cm.load_noise_files(
        glob.glob('/home/mkolopanis/psa64/21cmsense_noise/'
                  'dish_size_1/*drift_mod*.npz'))

noise_interp = py21cm.noise_interp2d(n_fqs, n_ks, noise)

kpls_pos = n.concatenate(n.array_split(kpls, [10, 11])[1:])
umag = 30/(3e8/(freq*1e9))
kpr = dk_du(z_bin) * umag
n_k = n.array(n.sqrt(kpls_pos**2 + kpr**2))
d2_n = noise_interp(freq*1e3, n_k)
pk_n = 2*n.pi**2/(n_k**3) * d2_n  # * 3887./2022
p_n = n.median(pk_n)
# p_n = n.max( pk_n )

pIvs_boot = []
pCvs_boot = []
pCs_boot = []
pIs_boot = []

print 'Boostrapping the Bootstraps'
for nboot in xrange(NBOOT):
    if (nboot+1) % 10 == 0:
        print '   ', nboot+1, '/', NBOOT
    dsum_pIvs = []
    dsum_pCvs = []
    dsum_pCs = []
    dsum_pIs = []
    # for cnt in xrange(int(n.ceil(pIvs.shape[-1]/96.))):
    for cnt in xrange(10):
        t = n.random.choice(range(pIvs.shape[-1]))
        b = n.random.choice(range(pIvs.shape[-2]))
        dsum_pIvs += [pIvs[:, :, b, t]]
        dsum_pCvs += [pCvs[:, :, b, t]]
        dsum_pIs += [pIs[:, :, b, t]]
        dsum_pCs += [pCs[:, :, b, t]]

    pIvs_sum = n.median(dsum_pIvs, 0)
    pCvs_sum = n.median(dsum_pCvs, 0)
    pIs_sum = n.median(dsum_pIs, 0)
    pCs_sum = n.median(dsum_pCs, 0)
    # Transpose arrays to make k's, injects
    pIvs_boot.append(n.ma.masked_invalid(pIvs_sum).T)
    pCvs_boot.append(n.ma.masked_invalid(pCvs_sum).T)
    pIs_boot.append(n.ma.masked_invalid(pIs_sum).T)
    pCs_boot.append(n.ma.masked_invalid(pCs_sum).T)
    # pIvs_neg, pIvs_0, pIvs_pos = n.array_split(pIvs_sum, [10,11],axis=-1)
    # pIvs_neg= pIvs_neg[::-1] #reorder in increaseing |kmag|

    # average pos and negative P(k) then concatenate P(0) with P(K_mag)
    # tmp = []
    # for cnt1 in xrange(2):
    #     h = n.random.choice(2)
    #     tmp = [pIvs_pos,pIvs_neg][h]
    # pIvs_fold= n.concatenate([pIvs_0,tmp],axis=1)
    # d2_boot.append( (pIvs_fold*n_k**3/(2*n.pi**2)).T )#makes kmag, inj
'''
# d2_boot = n.reshape(d2_boot,(n.shape(d2_boot)[1],
                      n.shape(d2_boot)[0]*n.shape(d2_boot)[-1]))
d2_pIv = n.mean(n.average(d2_boot,axis=0),-1)

d2_pIv_err = n.mean(errorbars(d2_boot,axis=0),-1)

d2_pIv_median = n.median(d2_boot,0).mean(-1)
d2_pIv_mean = n.mean(d2_boot,0).mean(-1)
d2_pIv_std = n.std(d2_boot,0).mean(-1)

p.figure()
# p.errorbar(n_k,d2_pIv_median,d2_pIv_err,fmt='kd',mfc='none',mec='black')
# p.errorbar(n_k,d2_pIv_mean,2*d2_pIv_std,fmt='rd',mfc='none',mec='red')
p.plot(n_k,n.abs(d2_pIv_median) + d2_pIv_err[1],'b--',label='pIv median')
p.plot(n_k,n.abs(d2_pIv_mean) + 2*d2_pIv_std,'g--',label='pIv mean')
p.plot(n_k,d2_n,'k--')
p.legend(loc='lower right')
p.grid(which='major')
p.show()
'''

pIvs_boot = n.ma.masked_invalid(pIvs_boot)
pCvs_boot = n.ma.masked_invalid(pCvs_boot)
pCs_boot = n.ma.masked_invalid(pCs_boot)
pIs_boot = n.ma.masked_invalid(pIs_boot)
# take the mean over bootstraps (which include eor draws and data bootstrap)
pIs = n.mean(pIs_boot, 0)
pCs = n.mean(pCs_boot, 0)
pCvs = n.mean(pCvs_boot, 0)
pCvs_pk = n.mean(pCvs_boot, 0)
pIvs_pk = n.mean(pIvs_boot, 0)
pIvs = n.mean(pIvs_boot, 0)

# find the 95% confidence interval over bootstraps
pI_err = errorbars(pIs_boot, axis=0)
pC_err = errorbars(pCs_boot, axis=0)
pCv_err = errorbars(pCvs_boot, axis=0)
pCv_pk_err = errorbars(pCvs_boot, axis=0)
pIv_pk_err = errorbars(pIvs_boot, axis=0)
pIv_err = errorbars(pIvs_boot, axis=0)

# all arrays should be of the form (num ks, num injects) now median over proper
# dimension
pIs_full = n.copy(pIs)
pCs_full = n.copy(pCs)
pCvs_full = n.copy(pCvs)
pIvs_full = n.copy(pIvs)

pI_err_full = n.copy(pI_err)
pC_err_full = n.copy(pC_err)
pCv_err_full = n.copy(pCv_err)
pIv_err_full = n.copy(pIv_err)


pIs = n.median(pIs, 0)
pCs = n.median(pCs, 0)
pCvs = n.median(pCvs, 0)
pIvs = n.median(pIvs, 0)
# median over injects to plot P(k)
pIvs_pk = n.median(pIvs_pk, -1)
pCvs_pk = n.median(pCvs_pk, -1)

# do same for errors
pI_err = n.median(pI_err, 1)
pC_err = n.median(pC_err, 1)
pCv_err = n.median(pCv_err, 1)
pIv_err = n.median(pIv_err, 1)

pIv_pk_err = n.median(pIv_pk_err, -1)
pCv_pk_err = n.median(pCv_pk_err, -1)


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
p.errorbar(pIs, pIvs, pIv_err, fmt='.', color='gray', alpha=.4)
p.errorbar(pIs, pCs, pC_err, fmt='k.')
p.plot(pIs, n.repeat(p_n, len(pIs)), 'g--')

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
"""p
for kpl,pk,err in zip(kpls,pks,errs):
    #p.([pklo,pkhi], [pk,pk], 'r')
    pkup = max(pk+err,1e-6)
    pkdn = max(pk-err,1e-6)
    print pkdn,pkup
    p.fill_between([pklo,pkhi], [pkdn,pkdn], [pkup,pkup], facecolor='gray',
                   edgecolor='gray')
"""

# Plot 3
ax3 = p.subplot(gs[5])  # used to be 3
p.setp(ax3.get_yticklabels(), visible=False)
# p.loglog(n.clip(pIs/pCs - 1, 1e-3,n.Inf), pCs, 'k.')
# p.loglog(n.abs(pIs/pCs - 1), pCs, 'k.')
p.errorbar((n.abs(pIs) + n.abs(p_n))/n.abs(pCs), n.abs(pCs),
           (n.abs(pIs) + n.abs(p_n))/(pCs)**2*pC_err +
           n.abs(pI_err)/n.abs(pCs),
           fmt='k.')
p.errorbar((n.abs(pIs) + n.abs(p_n))/n.abs(pCs), n.abs(pCvs),
           pCv_err, fmt='.', color='blue', alpha=.2)
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
"""
for kpl,pk,err in zip(kpls,pks,errs):
    pkup = max(pk+err,1e-6)
    pkdn = max(pk-err,1e-6)
    p.fill_between([1e-3,1e8], [pkdn,pkdn], [pkup,pkup], facecolor='gray',
                   edgecolor='gray')
    if pkup > 1e-6 and kpl > .17:
        sig_factors.append( sig_factor_interp(pkup))
"""

# Plot 1
ax0 = p.subplot(gs[1])  # used to be 0
p.setp(ax0.get_xticklabels(), visible=False)
p.errorbar(pIs, n.abs(pCs)/(n.abs(pIs) + n.abs(p_n)),
           pC_err/(n.abs(pIs) + n.abs(p_n)) + pCs/(n.abs(pIs) +
           n.abs(p_n))**2*pI_err, fmt='k.')
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
p.errorbar(kpls, pCvs_pk, pCv_pk_err, fmt='b.')
p.errorbar(kpls, pIvs_pk, pIv_pk_err, fmt='.', color='gray', alpha=.4)
ax0.set_xscale('log')
# p.errorbar(kpls, n.abs(pCvs), yerr=2*pCvs_err, fmt='k.', capsize=0)
p.grid()

p.savefig('sigloss_errorbar.png', format='png')

# Save p_out (lower, median, upper) and pCv (lower, median, upper) to file


labels = ['' for x in range(len(pCs)+1)]
labels[0] = 'P_data'
labels[1] = 'P_inj'
lowers = n.concatenate([[pCv_err[0][0]], pC_err[0]])
uppers = n.concatenate([[pCv_err[1][0]], pC_err[1]])
means = n.concatenate([[n.median(pCvs)], pCs])
_p_inj = n.concatenate([[0], pIs])
table = {'P_inj': _p_inj, 'Lower': lowers, 'Median': means, 'Upper': uppers}

ascii.write(table, 'sigloss_limits.dat',
            names=['P_inj', 'Lower', 'Median', 'Upper'])

pC_upper = pCs + pC_err[1]
pCv_upper = pCvs + pCv_err[1]

# p_eor = n.max([pIs[loc] for loc,val in enumerate(pC_tot)
#               if val <=pCv_tot[loc]])
p_inj_5 = interp1d(pC_upper, pIs, kind='linear', bounds_error=False,
                   fill_value=0, assume_sorted=False)
p_eor = p_inj_5(n.max(pCv_upper)).item()
beta = p_eor/n.max(pCv_upper)

sigma = n.median(pCv_err)/2
mu = n.median(pCvs)
ybins = n.logspace(n.log10(pklo), n.log10(pkhi), 500)
g_func = gauss(ybins, mu, sigma)
g_func /= g_func.sum()
g_interp = interp_gauss(ybins, g_func)
injects = (pCs + pC_err[1])
# inds = n.argsort(pIs)

pout_gauss = gauss(n.tile(ybins, (pIs_full.shape[0], pIs_full.shape[1], 1)).T,
                   pCs_full.T, n.median(pC_err_full, axis=0).T/2.)
pout_sum = pout_gauss.sum(0)
# pout_sum.shape += (1,)
pout_gauss /= pout_sum
dual_prob_func = pout_gauss.T * g_func
dual_prob = dual_prob_func.sum(-1)/(g_func**2).sum()

p_eor_lim = []

for cnt, pk in enumerate(dual_prob):
    inj_ind_95 = n.max(n.argwhere(n.log(pk) >= n.log10(.05)).squeeze())
    inj_ind_90 = n.max(n.argwhere(n.log(pk) >= n.log10(.1)).squeeze())

    p_eor_lim.append(pIs_full[cnt][inj_ind_95])
    # beta_95 = p_eor_95/pCs[inj_ind_95]
    # p_eor_90 = n.abs(pIs )[inj_ind_90]
    # beta_90 = p_eor_90/pCs[inj_ind_90]
'''
p.figure()
p.title('Total probabilty of $P_{inj}*P_{data}$/$P_{data} ^{2}$')
p.plot(pIs, dual_prob, 'k.')
p.xscale('log')
p.yscale('log')
p.xlabel('$P_{inj}$')
p.grid()
p.savefig('probability_plot.png',format='png')
p.figure()
p.title('Sample Probabilities $P_{inj}$')
p.plot(ybins, pout_gauss[:,10],'k.')
p.plot(ybins, pout_gauss[:,43],'g.')
p.plot(ybins, pout_gauss[:,49],'r.')
p.plot(ybins, g_func,'b.')
p.xscale('log')
p.yscale('log')
p.grid()

p.figure()
p.title('Sample Probabilities $P_{inj}*P_{data}$')
p.plot(ybins, dual_prob_func[10],'k.')
p.plot(ybins, dual_prob_func[43],'g.')
p.plot(ybins, dual_prob_func[49],'r.')
p.plot(ybins, g_func**2,'b.')
p.xscale('log')
p.yscale('log')
p.grid()
# prob = n.array([ g_interp(inj) for inj in injects ])

# beta = xb_fact_1[arg_050]
# beta_upper = xb_fact_1[arg_975] - beta
# beta_lower = beta - xb_fact_1[arg_025]
# #print ("Max sigloss factor z={0:.2f}:
#          {1:.2f}".format(z_bin,n.max(sig_factors)) )
f= open('sigloss_factor_confidence.txt','w')
f.write('95% confidence P_eor z={0:.2f}:  {1:.2e}\n'.format(z_bin,p_eor_95))
f.write('95% sigloss factor z={0:.2f}:  {1:.2f}\n'.format(z_bin,beta_95))
f.write('90% confidence P_eor z={0:.2f}:  {1:.2e}\n'.format(z_bin,p_eor_90))
f.write('90% sigloss factor z={0:.2f}:  {1:.2f}\n'.format(z_bin,beta_90))
f.close()
#
'''
with open('sigloss_p_eor_limits.txt', 'w') as f:
    f.write('k, Pk\n')
    for _k, _pk in zip(kpls, p_eor_lim):
        f.write('{0:.3e}, {1:.3e}\n'.format(_k, _pk))
f.close()

p.figure()
p.plot(kpls, p_eor_lim, 'r--')
p.errorbar(kpls, pCvs_pk, pCv_pk_err, fmt='k.')
p.grid()
p.yscale('log')
p.show()
