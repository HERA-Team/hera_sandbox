#! /usr/bin/env python
import numpy as n, pylab as p, sys
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import glob


### GETTING SIGLOSS DATA ###

pCs,pIs,pCvs = [],[],[]
pCs_full, pIs_full = [], []
pCs_err,pIs_err = [],[]
freq = []
for inject in glob.glob('inject_*'):
    print 'Reading', inject
    pspecs = glob.glob(inject + '/pspec_boot*.npz')
    inject = float(inject.split('_')[-1])
    pC_avg, pI_avg, pCv_avg = [], [], []
    pC_spec, pI_spec = [], []
    pcf,pif = [],[]
    for pspec in pspecs:
        npz = n.load(pspec)
        try: freq = npz['freq']
        except: pass
        pC,pI,pCv = npz['pk_vs_t'], npz['nocov_vs_t'], npz['pCv'] #(#chan, #times)
        kpls = n.array(npz['kpl'])
        pC_avg.append(n.average(pC.real)) #avg over freq and time
        pcf.append(pC.real)
        pI_avg.append(n.average(pI.real))
        pif.append(pI.real)
        #pCv_avg.append(n.average(pCv.real,axis=1)) #spectrum
        pCv_avg.append(pCv.real) #spectrum
        #pC_spec.append(n.average(pC.real, axis=1))
        #pI_spec.append(n.average(pI.real, axis=1))
        #pC_spec.append(pC.real)
        #pI_spec.append(pI.real)
    pCs.append(n.average(pC_avg)) #should bootstrap all these pC_avg's (20 of them) to get error bars... look at their distribution
    pCs_full.append(pcf)
    pIs_full.append(pif)
    pCs_err.append(n.std(pC_avg)/n.sqrt(len(pspecs)))
    pIs.append(n.average(pI_avg))
    pIs_err.append(n.std(pI_avg)/n.sqrt(len(pspecs)))
    #pCvs.append(n.average(pCv_avg,axis=0)) #spectrum
    pCvs.append(pCv_avg) #spectrum
    #pC_spec = n.average(pC_spec, axis=0)
    #print inject, pspec, pC_avg, pI_avg
    #p.figure(1)
    #p.subplot(211); p.loglog([pI_avg], [pC_avg], 'k.')
    #p.subplot(212); p.loglog([pC_avg], [pI_avg/pC_avg], 'k.')

pIs,pCs,pCvs = n.array(pIs), n.array(pCs), n.array(n.average(pCvs,axis=0)) #avg over inject #s
pIs_err,pCs_err = n.array(pIs_err), n.array(pCs_err)
pCs_full, pIs_full= n.array(pCs_full), n.array(pIs_full)
###Build an interpolator to find sigloss factors###
sig_factor_interp = interp1d(n.abs(pCs_full.ravel()), n.abs(pIs_full.ravel())/n.abs(pCs_full.ravel()),
                        kind='linear',bounds_error=False,fill_value=0)
order = n.argsort(n.abs(pCs)) #set up interpolator to work even if pCs are out of order
pCs_order = n.abs(pCs[order])
pIs_order = n.abs(pIs[order])
sig_factor_interp = interp1d(pCs_order, pIs_order/pCs_order,kind='linear',bounds_error=False,fill_value=0)

>>>>>>> carina/master
### GETTING PSPEC DATA ###
# XXX only used to get 'freq' variable

try:
    z_bin = f2z(freq)
except:
    print 'frequency not found in boots. Searching in pspec.npz'
    npz = n.load('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/lstbin_manybls/PS_frfnew/pspec_pk_k3pk.npz')
    #npz = n.load('/home/cacheng/capo/ctc/matt_data/noise_diffbls/pspec_pk_k3pk.npz') #matt's data
    freq = npz['freq']
    z_bin = f2z(freq)

#kpls,pks,errs = npz['kpl'], npz['pk'], npz['err']

### PLOTTING ###

fig = p.figure(1, figsize=(7,7))
gs = gridspec.GridSpec(2,3,width_ratios=[.4,1.2,.4],height_ratios=[.4,1]) #used to be 2,2
fig.subplots_adjust(left=.15, top=.95, bottom=.15, wspace=.35, hspace=.15, right=0.95)

#Plot 2
p.figure(1)
pklo,pkhi = 1e-4,1e10
ax2 = p.subplot(gs[4]) #used to be 2
#p.loglog(pIs, pCs, 'k.')
p.setp(ax2.get_yticklabels(), visible=False) #uncomment if no left-hand P(k) plot
p.errorbar(pIs, n.abs(pCs), xerr=2*pIs_err, yerr=2*pCs_err, capsize=0, fmt='k.')
p.loglog([pklo,pkhi],[pklo,pkhi], 'k-')
p.xlim(pklo, pkhi)
p.ylim(pklo, pkhi)
p.xlabel(r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize=14)
p.ylabel(r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize=14)
p.grid()
pkup = n.max(n.abs(pCvs))
pkdn = n.min(n.abs(pCvs))
p.fill_between([pklo,pkhi],[pkdn,pkdn],[pkup,pkup], facecolor='gray', edgecolor='gray')
"""
for kpl,pk,err in zip(kpls,pks,errs):
    #p.([pklo,pkhi], [pk,pk], 'r')
    pkup = max(pk+err,1e-6)
    pkdn = max(pk-err,1e-6)
    print pkdn,pkup
    p.fill_between([pklo,pkhi], [pkdn,pkdn], [pkup,pkup], facecolor='gray', edgecolor='gray')
"""

#Plot 3
ax3 = p.subplot(gs[5]) #used to be 3
p.setp(ax3.get_yticklabels(), visible=False)
#p.loglog(n.clip(pIs/pCs - 1, 1e-3,n.Inf), pCs, 'k.')
#p.loglog(n.abs(pIs/pCs - 1), pCs, 'k.')
p.errorbar(n.abs(pIs/pCs - 1), n.abs(pCs), xerr=2*pIs_err/pCs, yerr=2*pCs_err, fmt='k.', capsize=0)
print "pI/pC : ", pIs/pCs - 1
print "pI: ", pIs
print "pC: ", pCs
ax3.set_xscale('log')
ax3.set_yscale('log')
p.ylim(pklo, pkhi)
p.xlim(1e-1,1e5)
#p.xlim(1,10)
p.grid()
p.xlabel(r'$P_{\rm in}/P_{\rm out}-1$', fontsize=14)
p.fill_between([1e-3,1e8], [pkdn,pkdn], [pkup,pkup], facecolor='gray', edgecolor='gray')
"""
sig_factors = []
for kpl,pk,err in zip(kpls,pks,errs):
    pkup = max(pk+err,1e-6)
    pkdn = max(pk-err,1e-6)
    p.fill_between([1e-3,1e8], [pkdn,pkdn], [pkup,pkup], facecolor='gray', edgecolor='gray')
    if pkup > 1e-6 and kpl > .17:
        sig_factors.append( sig_factor_interp(pkup))
"""

#Plot 1
ax0 = p.subplot(gs[1]) #used to be 0
p.setp(ax0.get_xticklabels(), visible=False)
p.errorbar(pIs, n.abs(pCs/pIs), xerr=2*pIs_err, yerr=2*pCs_err/pIs, fmt='k.', capsize=0)
ax0.set_xscale('log')
p.xlim(pklo, pkhi)
#p.ylim(3e-2,3)
p.ylim(0,1.1)
p.grid()
p.ylabel(r'$P_{\rm out}/P_{\rm in}$', fontsize=14)

p.title('z = {0:.2f}'.format(z_bin))

#P(k) plot on left side
ax4 = p.subplot(gs[3])
p.setp(ax4.get_yticklabels(), visible=True)
ax4.set_xscale('log')
ax4.set_yscale('log')
p.ylim(pklo, pkhi)
p.plot(kpls,n.abs(pCvs.mean(0).mean(-1)),'k.')
#p.errorbar(kpls, n.abs(pCvs), yerr=2*pCvs_err, fmt='k.', capsize=0)
p.grid()

p.savefig('sigloss.png',format='png')

print sig_factors
print "Max sigloss factor z={0:.2f}:  {1:.2f}".format(z_bin,n.max(sig_factors))
f= open('sigloss_factor.txt','w')
f.write("Max sigloss factor z={0:.2f}:  {1:.2f}".format(z_bin,float(sig_factor_interp(n.max(n.abs(pCvs))))))
print "Max sigloss factor z={0:.2f}:  {1:.2f}".format(z_bin,float(sig_factor_interp(n.max(n.abs(pCvs))))) #n.max(sig_factors))
f.close()
p.show()



#p.show()
p.close()
