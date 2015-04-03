#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os, random
PLOT=True
NBOOT = 400
args = sys.argv[1:]

pks = {}
pk_2d = {}
temp_data = {}
nocov_2d = {}
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    kpl,cmd = f['kpl'], f['cmd']
    path = os.path.dirname(filename)
    if not pks.has_key(path):
        print '   ', path
        print '   ', cmd
        pks[path] = []
        pk_2d[path] = []
        temp_data[path] = []
        nocov_2d[path] = []
    pks[path].append(f['pk'])
    nancount =  n.sum(n.isnan(f['pk']))/float(f['pk'].size)
    if nancount>0:
        print 'Has NaN:',nancount,"skipping"
        continue
    pk_2d[path].append(f['pk_vs_t'])
    scalar = f['scalar']
    temp_data[path].append(f['temp_noise_var'])
    nocov_2d[path].append(f['nocov_vs_t'])
freq = f['freq']
paths = pks.keys()
print paths

pks = n.array([pks[path] for path in paths]) # (bltype,bootstraps,kpls)
pk_2d = n.array([pk_2d[path] for path in paths]) # (bltype,bootstraps,kpls,times)
temp_data = n.array([temp_data[path] for path in paths]) # (bltype,bootstraps,kpls,times), T averaged over all bls
nocov_2d = n.array([nocov_2d[path] for path in paths]) # (bltype,bootstraps,kpls,times), T averaged over all bls
nos_std_2d = n.var(temp_data, axis=1) * scalar # (bltype,kpls,times), thermal noise err in pk_2d
#tot_std_2d = n.std(pk_2d, axis=1) # (bltype,kpls,times), total err in pk_2d
avg_pk_2d = n.average(pk_2d, axis=1) # (bltype, kpls,times), best estimate of pk in each integration
#avg_pk_1d = n.average(avg_pk_2d, axis=2).real
wgts = n.ones_like(avg_pk_2d)
#avg_pk_1d.shape = avg_pk_1d.shape + (1,)
##wgts = 1./(nos_std_2d*tot_std_2d)
#wgts = 1./(nos_std_1d*(avg_pk_1d + nos_std_1d))
#nos_std_1d = n.median(nos_std_2d, axis=1)
#nos_std_1d.shape = (wgts.shape[0],1,wgts.shape[2])
#wgts = n.ones_like(wgts) / nos_std_1d**2
#wgts = 1./(n.abs(avg_pk_1d) + nos_std_2d)**2

if False: # override power spectrum with the version w/o covariance diagonalization
    print 'Overriding power spectrum with non-covariance diagonalized version'
    pk_2d = nocov_2d

#CLIP = True
CLIP = False
if CLIP:
    #pk_2d = pk_2d[...,250:550]
    #avg_pk_2d = avg_pk_2d[...,250:550]
    #wgts = wgts[...,250:550]
    pk_2d = pk_2d[...,200:800]
    avg_pk_2d = avg_pk_2d[...,200:800]
    wgts = wgts[...,200:800]
    ##pk_2d = pk_2d[...,200:600]
    ##wgts = wgts[...,200:600]
    #pk_2d = pk_2d[...,300:500]
    #wgts = wgts[...,300:500]
else:
    for i in xrange(nos_std_2d.shape[0]):
      for j in xrange(nos_std_2d.shape[1]):
        nos_std_2d[i,j] = n.convolve(nos_std_2d[i,j], n.ones((50,)), mode='same')
    wgts = 1./nos_std_2d**2

if PLOT: # plot some stuff
    plt1 = int(n.sqrt(len(paths)))
    plt2 = int(n.ceil(len(paths)/float(plt1)))
    for cnt,path in enumerate(paths):
        print cnt, path
        p.subplot(plt2,plt1,cnt+1)
        C.arp.waterfall(n.abs(n.average(temp_data[cnt], axis=0))**2 * scalar, mx=10, drng=3)
        p.colorbar(shrink=.5) 
    p.subplot(plt2,plt1,1); p.title(r'$|\langle\tilde V_b\rangle|^2$'); p.savefig('pspec_boot_1.png')
    p.clf()
    for cnt,path in enumerate(paths):
        p.subplot(plt2,plt1,cnt+1)
        C.arp.waterfall(nos_std_2d[cnt], mx=10, drng=3)
        p.colorbar(shrink=.5) 
    p.subplot(plt2,plt1,1); p.title('Thermal Noise [mK$^2$]'); p.savefig('pspec_boot_2.png')
    p.clf()
    for cnt,path in enumerate(paths):
        p.subplot(plt2,plt1,cnt+1)
        #C.arp.waterfall(avg_pk_2d[cnt], mx=10, drng=3)
        C.arp.waterfall(avg_pk_2d[cnt], mode='log', mx=8, drng=5)
        p.colorbar(shrink=.5) 
    p.subplot(plt2,plt1,1); p.title('Power Spectrum [mK$^2$]'); p.savefig('spec_boot_3.png')
    p.clf()
    plt1,plt2 = len(paths),3
    for cnt,path in enumerate(paths):
        p.subplot(plt2,plt1,0*plt1+cnt+1)
        #C.arp.waterfall(avg_pk_2d[cnt], mx=10, drng=4)
        C.arp.waterfall(avg_pk_2d[cnt], mode='real', mx=1e8, drng=2e8)
        p.colorbar(shrink=.5) 
        p.subplot(plt2,plt1,1*plt1+cnt+1)
        C.arp.waterfall(wgts[cnt])
        p.colorbar(shrink=.5) 
        p.subplot(plt2,plt1,2*plt1+cnt+1)
        #C.arp.waterfall(n.cumsum(avg_pk_2d[cnt]*wgts[cnt],axis=1)/n.cumsum(wgts[cnt],axis=1), mx=10, drng=4)
        C.arp.waterfall(n.cumsum(avg_pk_2d[cnt]*wgts[cnt],axis=1)/n.cumsum(wgts[cnt],axis=1), mode='real', mx=1e8, drng=2e8)
        p.colorbar(shrink=.5) 
    p.subplot(plt2,plt1,1); p.title('Weighted Power Spectrum [mK$^2$]'); p.savefig('pspec_boot_4.png')


    p.plot(avg_pk_2d[0,-1])
    p.plot(n.cumsum(avg_pk_2d[0,-1]*wgts[0,-1])/n.cumsum(wgts[0,-1]))
    p.savefig('pspec_boot_5.png')
print avg_pk_2d.shape, wgts.shape
print pk_2d.shape
print wgts.shape
pk_2d = pk_2d.transpose([1,2,3,0]).copy() # (bootstraps, kpls, times, bltypes)
pk_2d.shape = pk_2d.shape[:-2] + (pk_2d.shape[-2] * pk_2d.shape[-1],) # (bootstraps, kpls, timebls)
wgts = wgts.transpose([1,2,0]).copy() # (kpls, times, bltypes)
wgts.shape = wgts.shape[:-2] + (wgts.shape[-2] * wgts.shape[-1],) # (bootstraps, kpls, timebls)

#ntimes = pk_2d.shape[-1] / 2
ntimes = pk_2d.shape[-1]
print ntimes
pk_boot = []
pk_fold_boot = []
for boot in xrange(NBOOT):
    if boot % 10 == 0: print boot
    dsum,dwgt = 0, 0
    for t in xrange(ntimes):
        t = random.choice(range(pk_2d.shape[-1]))
        b = random.choice(range(pk_2d.shape[0]))
        #dsum += avg_pk_2d[:,t] * wgts[:,t]
        dsum += pk_2d[b,:,t] * wgts[:,t]
        dwgt += wgts[:,t]
    pk_boot.append(dsum/dwgt)
    k0 = n.abs(kpl).argmin()
    dsum_fold = dsum[k0:].copy()
    dwgt_fold = dwgt[k0:].copy()
    dsum_fold[1:] = 0
    dwgt_fold[1:] = 0
    dsum_pos,dwgt_pos = dsum[k0+1:].copy(), dwgt[k0+1:].copy()
    dsum_neg,dwgt_neg = dsum[k0-1:0:-1].copy(), dwgt[k0-1:0:-1].copy()
    for h in xrange(2): # bootstrap over which half of the spectrum (or both) are used
        h = random.randint(0,1)
        dsum_fold[1:] += [dsum_pos, dsum_neg][h]
        dwgt_fold[1:] += [dwgt_pos, dwgt_neg][h]
    pk_fold_boot.append(dsum_fold / dwgt_fold)
pk_boot = n.array(pk_boot).T
pk_fold_boot = n.array(pk_fold_boot).T

print 'Sorting bootstraps'
pk = n.average(pk_boot, axis=1)
pk_fold = n.average(pk_fold_boot, axis=1)
# this is excluding imag component in noise estimate `
pk_boot = n.sort(pk_boot.real, axis=1) # losing imag component here
pk_fold_boot = n.sort(pk_fold_boot.real, axis=1) # losing imag component here
#pk_boot = n.sort(pk_boot.imag, axis=1) # losing imag component here
#pk_fold_boot = n.sort(pk_fold_boot.imag, axis=1) # losing imag component here
if True:
    print 'Deriving errors from histogram'
    up_thresh = int(n.around(0.975 * pk_boot.shape[1])) # 2 sigma, single tail
    # important to only include real component in estimation of error
    err = (pk_boot[:,up_thresh] - pk.real) / 2 # effective "1 sigma" derived from actual 2 sigma
#    err = (pk_boot[:,up_thresh] - pk.imag) / 2 # effective "1 sigma" derived from actual 2 sigma
    up_thresh_fold = int(n.around(0.975 * pk_fold_boot.shape[1])) # 2 sigma, single tail
    err_fold = (pk_fold_boot[:,up_thresh_fold] - pk_fold.real) / 2 # effective "1 sigma" derived from actual 2 sigma
#    err_fold = (pk_fold_boot[:,up_thresh_fold] - pk_fold.imag) / 2 # effective "1 sigma" derived from actual 2 sigma
    for k_,pk_,err_ in zip(kpl, pk/1e6, 2*err/1e6):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (k_, pk_+err_,pk_,err_)
    print '-'*70
    for k_,pk_,err_ in zip(kpl[k0:], pk_fold/1e6, 2*err_fold/1e6):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (k_, pk_+err_,pk_,err_)
else:
    err = n.std(pk_boot, axis=1)
    err_fold = n.std(pk_fold_boot, axis=1)
if PLOT:
    #p.subplot(221); C.arp.waterfall(pk_boot, mx=9, drng=3); p.colorbar(shrink=.5)
    p.subplot(221); p.plot(pk_boot)
    p.subplot(222)
    p.plot(pk)
    p.plot(n.median(pk_boot,axis=1))
    p.plot(pk+2*err)
    #p.subplot(223); C.arp.waterfall(pk_fold_boot, mx=9, drng=3); p.colorbar(shrink=.5)
    p.subplot(223); p.plot(pk_fold_boot)
    p.subplot(224)
    p.plot(pk_fold)
    p.plot(n.median(pk_fold_boot,axis=1))
    p.plot(pk_fold+2*err_fold)
    p.savefig('pspec_boot_6.png')

print 'Writing pspec.npz'
print "saving freq=",freq
n.savez('pspec.npz', kpl=kpl, pk=pk, err=err, pk_fold=pk_fold, err_fold=err_fold, cmd=cmd,freq=freq)
    
    

#if True: # new way, redo time-dependent weighting
#    #pk_2d = n.array(pk_2d).real
#    pks = n.array(pks)
#    pk_2d = n.array(pk_2d)
#    avg_2d = n.average(pk_2d, axis=0)
#    nos_var_2d = n.var(temp_noise_var, axis=0) * scalar
#    #temp_noise_var_2d = n.average(temp_noise_var_2d, axis=0)
#    #temp_noise_var_2d.shape = (1,) + temp_noise_var_2d.shape
#    std_2d = n.std(pk_2d, axis=0)
#    var_2d = n.var(pk_2d, axis=0)
#    sig_var_2d = 0.5*(var_2d/nos_var_2d - nos_var_2d)
#    #wgt_2d = 1/std_2d**2
#    wgt_2d_a = 1/var_2d
#    wgt_2d_b = 1/nos_var_2d**2
#    #wgt_2d_c = sig_var_2d**2/(sig_var_2d**2 + var_2d)
#    wgt_2d_c = avg_2d**2/(avg_2d**2 + var_2d)
#    #wgt_2d_c = 1/n.median(nos_var_2d**2, axis=0)
#    #wgt_2d_c = n.convolve(n.ones(80), wgt_2d_c, mode='same')
#    #wgt_2d_c.shape = (1,) + wgt_2d_c.shape
#
#    p.subplot(3,4, 1); C.arp.waterfall(avg_2d,               drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4, 2); C.arp.waterfall(var_2d,        mx=19, drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4, 3); C.arp.waterfall(nos_var_2d**2, mx=19, drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4, 4); C.arp.waterfall(sig_var_2d**2, mx=19, drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4, 6); C.arp.waterfall(wgt_2d_a,    drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4, 7); C.arp.waterfall(wgt_2d_b,    drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4, 8); C.arp.waterfall(wgt_2d_c,    drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4,10); C.arp.waterfall(avg_2d * wgt_2d_a,    drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4,11); C.arp.waterfall(avg_2d * wgt_2d_b,    drng=4); p.colorbar(shrink=.5)
#    p.subplot(3,4,12); C.arp.waterfall(avg_2d * wgt_2d_c,    drng=4); p.colorbar(shrink=.5)
#    p.show()
#    
#    wgt_2d_a.shape = (1,) + wgt_2d_a.shape
#    wgt_2d_b.shape = (1,) + wgt_2d_b.shape
#    wgt_2d_c.shape = (1,) + wgt_2d_c.shape
#    pks = n.average(pk_2d, axis=2)
#    pks_a = n.sum(pk_2d * wgt_2d_a, axis=2) / n.sum(wgt_2d_a, axis=2)
#    pks_b = n.sum(pk_2d * wgt_2d_b, axis=2) / n.sum(wgt_2d_b, axis=2)
#    pks_c = n.sum(pk_2d * wgt_2d_c, axis=2) / n.sum(wgt_2d_c, axis=2)
#
##combine = n.average
#combine = n.median
#pk = combine(pks.real, axis=0)
#err = n.std(pks.real, axis=0)
##C.arp.waterfall(pks.real); p.colorbar(shrink=.5); p.show()
#pk_a = combine(pks_a.real, axis=0)
#err_a = n.std(pks_a.real, axis=0)
#pk_b = combine(pks_b.real, axis=0)
#err_b = n.std(pks_b.real, axis=0)
#pk_c = combine(pks_c.real, axis=0)
#err_c = n.std(pks_c.real, axis=0)
#
#p.plot(kpl, pk, 'k.-')
#p.plot(kpl, pk_a, 'c.-')
#p.plot(kpl, pk_b, 'm.-')
#p.plot(kpl, pk_c, 'g.-')
#p.show()
##err = n.sqrt(n.median(n.abs(pks-pk)**2, axis=0))
#
#pk,err = pk_c, err_c

#print 'Writing pspec.npz'
#n.savez('pspec.npz', kpl=kpl, pk=pk, err=err, cmd=cmd)

