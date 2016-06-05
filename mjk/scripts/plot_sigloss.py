# ! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import numpy as n, pylab as p
from matplotlib import gridspec
import glob
import sys, optparse
import ipdb
from capo.pspec import f2z
from scipy.interpolate import interp1d

fig = p.figure(1, figsize=(7,7))
gs = gridspec.GridSpec(2,2, width_ratios=[1,.4], height_ratios=[.4,1])
fig.subplots_adjust(left=.15, top=.95, bottom=.15, wspace=.15, hspace=.15, right=0.95)

o=optparse.OptionParser()
o.add_option('--path',dest='path', action='store', type='string',\
        default='',help='Provide output folder to find power spectrum for comparison')
o.add_option('--pspec',dest='pspec', action='store', type='string',\
        default='',help='Provide output prefix to find power spectrum for comparison')
o.add_option('--pol',dest='pol', action='store', type='string',\
        default='',help='Provide polarization prefix to find power spectrum for comparison')
opts,args=o.parse_args(sys.argv[1:])
#fig2 = p.figure(2)
pCs,pIs = [],[]
pCs_err,pIs_err = [],[]
chan=[]
for inject in args:
    chan = inject.split('/')[1]
    pspecs = glob.glob(inject + '/pspec_boot*.npz')
    pspecs.sort()
    #ipdb.set_trace()
    inject = float(inject.split('_')[-1])
    pC_avg, pI_avg = [], []
    pC_spec, pI_spec = [], []
    for pspec in pspecs:
        npz = n.load(pspec)
        pC,pI = npz['pk_vs_t'], npz['nocov_vs_t']
        pC_avg.append(n.average(pC.real))
        pI_avg.append(n.average(pI.real))
        #pC_spec.append(n.average(pC.real, axis=1))
        #pI_spec.append(n.average(pI.real, axis=1))
        #pC_spec.append(pC.real)
        #pI_spec.append(pI.real)
    pCs.append(n.average(pC_avg))
    pCs_err.append(n.std(pC_avg)/n.sqrt(len(pspecs)))
    pIs.append(n.average(pI_avg))
    pIs_err.append(n.std(pI_avg)/n.sqrt(len(pspecs)))
    #pC_spec = n.average(pC_spec, axis=0)
    #pI_spec = n.average(pI_spec, axis=0)
    #p.figure(2); p.plot(pI_spec/pC_spec)
    #print inject, pspec, pC_avg, pI_avg
    #p.figure(1)
    #p.subplot(211); p.loglog([pI_avg], [pC_avg], 'k.')
    #p.subplot(212); p.loglog([pC_avg], [pI_avg/pC_avg], 'k.')
#p.show()
pIs,pCs = n.array(pIs), n.array(pCs)
pIs_err,pCs_err = n.array(pIs_err), n.array(pCs_err)
p.figure(1)

pklo,pkhi = 1e-3,1e9
p.subplot(gs[2])
#p.loglog(pIs, pCs, 'k.')
p.errorbar(pIs, pCs, xerr=2*pIs_err, yerr=2*pCs_err, capsize=0, fmt='k.')
p.loglog([pklo,pkhi],[pklo,pkhi], 'k-')
p.xlim(pklo, pkhi)
p.ylim(pklo, pkhi)
p.xlabel(r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize=14)
p.ylabel(r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize=14)
p.grid()


#npz = n.load('/home/mkolopanis/psa64/psa64_multiz/data/pspec_Nov18_vanilla_skew_frpad_2_'+chan+'_I.npz')

sig_factor_interp = interp1d(pIs, pIs/pCs,kind='linear',bounds_error=False,fill_value=0)

npz_file = opts.path+'/pspec_'+opts.pspec+'_'+chan+'_'+opts.pol+'.npz'
npz = n.load(npz_file)
#npz = n.load('/home/mkolopanis/psa64/pspec_Mar17_vanilla_ali_frf_recon_'+chan+'_I.npz')
#print npz.files
kpls,pks,errs = npz['kpl'], npz['pk'], npz['err']
freq = npz['freq']
z_bin = f2z(freq)

sig_factors= []
for kpl,pk,err in zip(kpls,pks,errs):
    #p.([pklo,pkhi], [pk,pk], 'r')
    pkup = max(pk+err,1e-6)
    pkdn = max(pk-err,1e-6)
    p.subplot(gs[2])
    p.fill_between([pklo,pkhi], [pkdn,pkdn], [pkup,pkup], facecolor='gray', edgecolor='gray')
    p.subplot(gs[3])
    p.fill_between([1e-3,20], [pkdn,pkdn], [pkup,pkup], facecolor='gray', edgecolor='gray')
    if pkup > 1e-6:
        sig_factors.append( sig_factor_interp(pkup))

sig_factors = n.array(sig_factors)

ax3 = p.subplot(gs[3])
p.setp(ax3.get_yticklabels(), visible=False)
#p.loglog(n.clip(pIs/pCs - 1, 1e-3,n.Inf), pCs, 'k.')
#p.loglog(n.abs(pIs/pCs - 1), pCs, 'k.')
p.errorbar(pIs/pCs - 1, pCs, xerr=2*pIs_err/pCs, yerr=2*pCs_err, fmt='k.', capsize=0)
#print pIs/pCs - 1
#print pCs
ax3.set_xscale('log')
ax3.set_yscale('log')
p.ylim(pklo, pkhi)
p.xlim(1e-3,20)
#p.xlim(1,10)
p.grid()
p.xlabel(r'$P_{\rm in}/P_{\rm out}-1$', fontsize=14)

ax0 = p.subplot(gs[0])
p.setp(ax0.get_xticklabels(), visible=False)
p.errorbar(pIs, pCs/pIs, xerr=2*pIs_err, yerr=2*pCs_err/pIs, fmt='k.', capsize=0)
ax0.set_xscale('log')
p.xlim(pklo, pkhi)
#p.ylim(3e-2,3)
p.ylim(0,1.1)
p.grid()
p.ylabel(r'$P_{\rm out}/P_{\rm in}$', fontsize=14)


p.title('z = {0:.2f}'.format(z_bin))

p.savefig('sigloss.png',format='png')

print "Max sigloss factor z={0:.2f}:  {1:.2f}".format(z_bin,n.max(sig_factors))

p.show()


