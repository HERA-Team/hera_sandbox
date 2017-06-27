#! /usr/bin/env python
import numpy as n, pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import glob
import optparse, sys

o = optparse.OptionParser()
o.add_option('--output', type='string', default=None,
    help='Output filename to save power spectrum values.')
opts,args = o.parse_args(sys.argv[1:])


### GETTING SIGLOSS DATA ###

pCs,pIs,pCvs,pIvs = [],[],[],[]
pCs_err,pIs_err = [],[]
freq = []
for inject in glob.glob('inject_*'):
    print 'Reading', inject
    pspecs = glob.glob(inject + '/pspec_boot*.npz') 
    inject = float(inject.split('_')[-1])
    pC_avg, pI_avg, pCv_avg, pIv_avg = [], [], [], []
    pC_spec, pI_spec = [], []
    for pspec in pspecs:
        npz = n.load(pspec)
        try: freq = npz['freq']
        except: pass
        pC,pI = npz['pk_vs_t'], npz['nocov_vs_t']
        try: pCv,pIv =  npz['pCv'], npz['pIv'] #(#chan, #times)
        except: pass
        kpls = n.array(npz['kpl'])
        pC_avg.append(n.average(pC.real)) #avg over freq and time
        pI_avg.append(n.average(pI.real))
        try: 
            pCv_avg.append(n.average(pCv.real,axis=1)) #spectrum
            pIv_avg.append(n.average(pIv.real,axis=1))
        except: pass
        #pC_spec.append(n.average(pC.real, axis=1))
        #pI_spec.append(n.average(pI.real, axis=1))
        #pC_spec.append(pC.real)
        #pI_spec.append(pI.real)
    pCs.append(n.average(pC_avg)) 
    pCs_err.append(n.std(pC_avg)/n.sqrt(len(pspecs)))
    pIs.append(n.average(pI_avg))
    pIs_err.append(n.std(pI_avg)/n.sqrt(len(pspecs)))
    try:  
        pCvs.append(n.average(pCv_avg,axis=0)) #spectrum
        pIvs.append(n.average(pIv_avg,axis=0))
    except: pass
    #pC_spec = n.average(pC_spec, axis=0)
    #print inject, pspec, pC_avg, pI_avg
    #p.figure(1)
    #p.subplot(211); p.loglog([pI_avg], [pC_avg], 'k.')
    #p.subplot(212); p.loglog([pC_avg], [pI_avg/pC_avg], 'k.')

pIs,pCs,pCvs,pIvs = n.array(pIs), n.array(pCs), n.array(n.average(pCvs,axis=0)), n.array(n.average(pIvs,axis=0)) #avg over inject #s
pIs_err,pCs_err = n.array(pIs_err), n.array(pCs_err)

#Plot pCv points from pspec_cov_v004 output, not sigloss output
pCv_points = n.load('pspec_C.npz')['pk']
pCvs = n.abs(pCv_points) #XXX overwrites pCvs from inject outputs
pIv_points = n.load('pspec_I.npz')['pk']
pIvs = n.abs(pIv_points)

###Build an interpolator to find sigloss factors###
order = n.argsort(n.abs(pCs)) #set up interpolator to work even if pCs are out of order
pCs_order = n.abs(pCs[order])
pIs_order = n.abs(pIs[order])
try: sig_factor_interp = interp1d(pCs_order, pIs_order/pCs_order,kind='linear',bounds_error=False,fill_value=0)
except: pass

#####  NOTE: the following came from capo/arp in conflict with an empty line. I have left it here but commented.
#<BEGIN COMMENT>
#=======
#### GETTING PSPEC DATA ###
## XXX only used to get 'freq' variable
#
#try:
#    z_bin = f2z(freq)
#except:
#    print 'frequency not found in boots. Searching in pspec.npz'
#    npz = n.load('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/lstbin_manybls/PS_frfnew/pspec_pk_k3pk.npz')
#    #npz = n.load('/home/cacheng/capo/ctc/matt_data/noise_diffbls/pspec_pk_k3pk.npz') #matt's data
#    freq = npz['freq']
#    z_bin = f2z(freq)
#
##kpls,pks,errs = npz['kpl'], npz['pk'], npz['err']
#
#### PLOTTING ###
#
#fig = p.figure(1, figsize=(7,7))
#gs = gridspec.GridSpec(2,3,width_ratios=[.4,1.2,.4],height_ratios=[.4,1]) #used to be 2,2
#fig.subplots_adjust(left=.15, top=.95, bottom=.15, wspace=.35, hspace=.15, right=0.95)
#
##Plot 2
#p.figure(1)
#pklo,pkhi = 1e-6,1e18 #1e2,1e14
#ax2 = p.subplot(gs[4]) #used to be 2
##p.loglog(pIs, pCs, 'k.')
#p.setp(ax2.get_yticklabels(), visible=False) #uncomment if no left-hand P(k) plot
#p.errorbar(pIs, n.abs(pCs), xerr=2*pIs_err, yerr=2*pCs_err, capsize=0, fmt='k.')
#p.loglog([pklo,pkhi],[pklo,pkhi], 'k-')
#p.xlim(pklo, pkhi)
#p.ylim(pklo, pkhi)
#p.xlabel(r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize=14)
#p.ylabel(r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize=14)
#p.grid()
#pkup = max(n.abs(pCvs))
#pkdn = min(n.abs(pCvs))
#p.fill_between([pklo,pkhi],[pkdn,pkdn],[pkup,pkup], facecolor='gray', edgecolor='gray')
#"""
#for kpl,pk,err in zip(kpls,pks,errs):
#    #p.([pklo,pkhi], [pk,pk], 'r')
#    pkup = max(pk+err,1e-6)
#    pkdn = max(pk-err,1e-6)
#    print pkdn,pkup
#    p.fill_between([pklo,pkhi], [pkdn,pkdn], [pkup,pkup], facecolor='gray', edgecolor='gray')
#"""
#
##Plot 3    
#ax3 = p.subplot(gs[5]) #used to be 3
#p.setp(ax3.get_yticklabels(), visible=False)
##p.loglog(n.clip(pIs/pCs - 1, 1e-3,n.Inf), pCs, 'k.')
##p.loglog(n.abs(pIs/pCs - 1), pCs, 'k.')
#p.errorbar(n.abs(pIs/pCs - 1), n.abs(pCs), xerr=2*pIs_err/pCs, yerr=2*pCs_err, fmt='k.', capsize=0)
#print "pI/pC : ", pIs/pCs - 1
#print "pI: ", pIs
#print "pC: ", pCs
#ax3.set_xscale('log')
#ax3.set_yscale('log')
#p.ylim(pklo, pkhi)
#p.xlim(1e-1,1e5)
##p.xlim(1,10)
#p.grid()
#p.xlabel(r'$P_{\rm in}/P_{\rm out}-1$', fontsize=14)
#p.fill_between([1e-3,1e8], [pkdn,pkdn], [pkup,pkup], facecolor='gray', edgecolor='gray')
#"""
#sig_factors = []
#sig_factors.append(sig_factor_interp(pkup))
#for kpl,pk,err in zip(kpls,pks,errs):
#    pkup = max(pk+err,1e-6)
#    pkdn = max(pk-err,1e-6)
#    p.fill_between([1e-3,1e8], [pkdn,pkdn], [pkup,pkup], facecolor='gray', edgecolor='gray')
#    if pkup > 1e-6 and kpl > .17:
#        sig_factors.append( sig_factor_interp(pkup))
#"""
#
##Plot 1
#ax0 = p.subplot(gs[1]) #used to be 0
#p.setp(ax0.get_xticklabels(), visible=False)
#p.errorbar(pIs, n.abs(pCs/pIs), xerr=2*pIs_err, yerr=2*pCs_err/pIs, fmt='k.', capsize=0)
#ax0.set_xscale('log')
#p.xlim(pklo, pkhi)
##p.ylim(3e-2,3)
#p.ylim(0,1.1)
#p.grid()
#p.ylabel(r'$P_{\rm out}/P_{\rm in}$', fontsize=14)
#
#p.title('z = {0:.2f}'.format(z_bin))
#p.savefig('sigloss.png',format='png')
#
##P(k) plot on left side
#ax4 = p.subplot(gs[3])
#p.setp(ax4.get_yticklabels(), visible=True)
#ax4.set_xscale('log')
#ax4.set_yscale('log')
#p.ylim(pklo, pkhi)
#<END COMMENT>
p.plot(n.abs(kpls),n.abs(pCvs),'k.')
#p.errorbar(kpls, n.abs(pCvs), yerr=2*pCvs_err, fmt='k.', capsize=0)
p.grid()

try:
    for k,kpl in enumerate(kpls):
        print ("%5.5f" % kpl), ':', ("%5.5f" % n.abs(pCvs[k])), ':', float(sig_factor_interp(n.abs(pCvs[k])))
    maxfactor = float(sig_factor_interp(n.max(n.abs(pCvs))))
    print "Max sigloss factor z={0:.2f}:  {1:.2f}".format(z_bin,maxfactor) #n.max(sig_factors))
except: pass

p.show()

if opts.output != None:
    print 'Writing '+opts.output+'.npz'
    n.savez(opts.output+'.npz', kpl=kpls, pCv=pCvs, pIv=pIvs, factor=maxfactor, cmd=' '.join(sys.argv))

