#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os, random
#sample run: python pspec_boot_v2.py boot-11/* boot01/* boot11/* coot-1101/* coot1101/*
o = optparse.OptionParser()
o.add_option('--nocov', action='store_true', 
            help='Use non covariance applied data')
o.add_option('--odir',default='out/', action='store',
            help='name of output directory')
o.add_option('--fname',default='pspec', action='store',
            help='Name of output file')
opts,args = o.parse_args(sys.argv[1:])

OutDIR = opts.odir
fname = opts.fname
NBOOT = 400
MEDIAN = True
CLIP = False
LO,HI = 40,320
NOCOV=opts.nocov
#LO,HI = 40,600

pk_vs_t = {}
err_vs_t = {}
temp_noise_dic = {}
nocov_vs_t = {}
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    kpl,cmd = f['kpl'], f['cmd']
    path = os.path.dirname(filename)
    if not pk_vs_t.has_key(path):
        print '   ', path
        print '   ', cmd
        pk_vs_t[path] = []
        err_vs_t[path] = []
        temp_noise_dic[path] = []
        nocov_vs_t[path] = []
    pk_vs_t[path].append(f['pk_vs_t'])
    scalar = f['scalar']
    #err_vs_t[path].append(n.average(f['err_vs_t'][:,120:141], axis=1))
    #temp_noise_var[path].append(n.average(f['temp_noise_var'][:,120:141], axis=1))
    nocov_vs_t[path].append(f['nocov_vs_t'])

paths = pk_vs_t.keys()
k0 = n.abs(kpl).argmin()

##############################################################################
#NN = n.min([len(pk_vs_t['boot'][i][0]) for i in range(len(pk_vs_t['boot']))])
#print NN

#for i in range(len(pk_vs_t['boot'])):
#    for j in range(len(pk_vs_t['boot'][i])):
#        pk_vs_t['boot'][i][j] = pk_vs_t['boot'][i][j][:NN]

pk_2d = n.empty(len(paths), dtype=n.ndarray)
nocov_2d = n.empty(len(paths), dtype=n.ndarray)
temp_noise_var = n.empty(len(paths), dtype=n.ndarray)
for i,path in enumerate(paths):
    pk_2d[i] = n.array(pk_vs_t[path]) # (bltype/path,bootstraps/numfilesinpath,kpls/nchan,times)
    nocov_2d[i] = n.array(nocov_vs_t[path]) # (bltype,bootstraps,kpls,times), T averaged over all bls
    temp_noise_var[i] = n.array(temp_noise_dic[path])
#err_2d = n.array([err_vs_t[path] for path in paths])
#print err_2d.shape, temp_noise_var.shape, pk_2d.shape
#wgts = 1./(temp_noise_var * err_2d)
#wgts.shape = wgts.shape[:2] + (1,) + wgts.shape[2:]
#import IPython; IPython.embed()
#average over files in each path
avg_pk_2d = [n.average(pk_2d[i], axis=0) for i, k in enumerate(pk_2d)]
wgts = [n.ones_like(avg_pk_2d[i]) for i,av in enumerate(avg_pk_2d)]
#avg_pk_2d, wgts = n.array([]),n.array([])
#for i, pk in enumerate(pk_2d):
#    n.append(avg_pk_2d,n.average(pk_2d[i], axis=0))# (bltype, kpls,times), best estimate of pk in each integration
#    n.append(wgts,n.ones_like(n.average(pk_2d[i], axis=0)))
#import IPython; IPython.embed()
#wgts = n.ones_like(avg_pk_2d)
#avg_pk_1d.shape = avg_pk_1d.shape + (1,)
##wgts = 1./(nos_std_2d*tot_std_2d)
#wgts = 1./(nos_std_1d*(avg_pk_1d + nos_std_1d))
#nos_std_1d = n.median(nos_std_2d, axis=1)
#nos_std_1d.shape = (wgts.shape[0],1,wgts.shape[2])
#wgts = n.ones_like(wgts) / nos_std_1d**2
#wgts = 1./(n.abs(avg_pk_1d) + nos_std_2d)**2

if NOCOV: # override power spectrum with the version w/o covariance diagonalization
    print 'Overriding power spectrum with non-covariance diagonalized version'
    pk_2d = nocov_2d

if False: # XXX decimate
    DEC = 4
    pk_2d = pk_2d[...,::DEC]
    wgts = wgts[...,::DEC]
    avg_pk_2d = avg_pk_2d[...,::DEC]

if CLIP:
    pk_2d = pk_2d[...,LO:HI]
    avg_pk_2d = avg_pk_2d[...,LO:HI]
    wgts = wgts[...,LO:HI]
else:
    pass
    #for i in xrange(nos_std_2d.shape[0]):
    #  for j in xrange(nos_std_2d.shape[1]):
    #    nos_std_2d[i,j] = n.convolve(nos_std_2d[i,j], n.ones((50,)), mode='same')
    #wgts = 1./nos_std_2d**2

#if True: # plot some stuff
if False: # plot some stuff
    plt1 = int(n.sqrt(len(paths)))
    plt2 = int(n.ceil(len(paths)/float(plt1)))
    #for cnt,path in enumerate(paths):
    #    print cnt, path
    #    p.subplot(plt2,plt1,cnt+1)
    #    C.arp.waterfall(n.abs(n.average(temp_data[cnt], axis=0))**2 * scalar, mx=10, drng=3)
    #    p.colorbar(shrink=.5) 
    #p.subplot(plt2,plt1,1); p.title(r'$|\langle\tilde V_b\rangle|^2$'); p.show()
    #for cnt,path in enumerate(paths):
    #    p.subplot(plt2,plt1,cnt+1)
    #    C.arp.waterfall(nos_std_2d[cnt], mx=10, drng=3)
    #    p.colorbar(shrink=.5) 
    #p.subplot(plt2,plt1,1); p.title('Thermal Noise [mK$^2$]'); p.show()
    for cnt,path in enumerate(paths):
        p.subplot(plt2,plt1,cnt+1)
        #C.arp.waterfall(avg_pk_2d[cnt], mx=10, drng=3)
        C.arp.waterfall(avg_pk_2d[cnt], mode='real', mx=5e7, drng=1e8)
        p.colorbar(shrink=.5) 
    p.subplot(plt2,plt1,1); p.title('Power Spectrum [mK$^2$]'); p.show()
    plt1,plt2 = len(paths),3
    for cnt,path in enumerate(paths):
        p.subplot(plt2,plt1,0*plt1+cnt+1)
        #C.arp.waterfall(avg_pk_2d[cnt], mx=10, drng=4)
        C.arp.waterfall(avg_pk_2d[cnt], mode='real', mx=5e7, drng=1e8)
        p.colorbar(shrink=.5) 
        p.subplot(plt2,plt1,1*plt1+cnt+1)
        C.arp.waterfall(wgts[cnt])
        p.colorbar(shrink=.5) 
        p.subplot(plt2,plt1,2*plt1+cnt+1)
        #C.arp.waterfall(n.cumsum(avg_pk_2d[cnt]*wgts[cnt],axis=1)/n.cumsum(wgts[cnt],axis=1), mx=10, drng=4)
        C.arp.waterfall(n.cumsum(avg_pk_2d[cnt]*wgts[cnt],axis=1)/n.cumsum(wgts[cnt],axis=1), mode='real', mx=5e7, drng=1e8)
        p.colorbar(shrink=.5) 
    p.subplot(plt2,plt1,1); p.title('Weighted Power Spectrum [mK$^2$]'); p.show()

#print avg_pk_2d.shape, wgts.shape
for i in range(len(wgts)): print 'avg_pk_2d[i].shape, wgts[i].shape',avg_pk_2d[i].shape, wgts[i].shape

#p.plot(avg_pk_2d[0,30])
#p.plot(n.cumsum(avg_pk_2d[0,30]*wgts[0,30])/n.cumsum(wgts[0,30]))
#p.show()

#print pk_2d.shape
#print wgts.size
#pk_2d = pk_2d.transpose([1,2,3,0]).copy() # (bootstraps(from pspec_cov4), kpls, times, bltypes(number paths))
#pk_2d.shape = pk_2d.shape[:-2] + (pk_2d.shape[-2] * pk_2d.shape[-1],) # (bootstraps, kpls, time*bltypes)
#gts = wgts.transpose([1,2,0]).copy() # (kpls, times, bltypes)
#gts.shape = wgts.shape[:-2] + (wgts.shape[-2] * wgts.shape[-1],) # (bootstraps, kpls, time*bltypes)
pk_2d = n.concatenate(tuple(pk_2d),axis=-1)
wgts = n.concatenate(tuple(wgts),axis=-1)
print 'pk_2d.shape',pk_2d.shape
print 'wgts.shape',wgts.shape
########################################################################################################
#ntimes = pk_2d.shape[-1] / 2
npaths = pk_2d.shape[0]
ntimes = pk_2d.shape[-1]  #now its time*bls
print 'ntimes, npaths=  ', ntimes,npaths
pk_boot = []
pk_fold_boot = []
for boot in xrange(NBOOT):
    if boot % 10 == 0: print boot
    #dsum,dwgt = 0, 0
    dsum,dwgt = [],[]
    dsum_fold,dwgt_fold = [], []
    for t in xrange(ntimes):
        t = random.choice(range(pk_2d.shape[-1]))
        b = random.choice(range(pk_2d.shape[0]))
        #dsum += pk_2d[b,:,t] * wgts[:,t]
        #dwgt += wgts[:,t]
        dsum += [pk_2d[b,:,t] * wgts[:,t]]
        dwgt += [wgts[:,t]]
    for t in xrange(2*ntimes):
        t = random.choice(range(pk_2d.shape[-1]))
        b = random.choice(range(pk_2d.shape[0]))
        h = random.randint(0,1)
        if h == 0:
            dsum_fold += [pk_2d[b,k0+1:,t] * wgts[k0+1:,t]]
            dwgt_fold += [wgts[k0+1:,t]]
        else:
            dsum_fold += [pk_2d[b,k0-1::-1,t] * wgts[k0-1::-1,t]]
            dwgt_fold += [wgts[k0-1::-1,t]]
    if MEDIAN:
        dsum,dwgt = n.median(dsum, axis=0), n.median(dwgt, axis=0)
        dsum_fold,dwgt_fold = n.median(dsum_fold, axis=0), n.median(dwgt_fold, axis=0)
    else:
        dsum,dwgt = n.average(dsum,axis=0), n.average(dwgt,axis=0)
    #_dsum_fold,_dwgt_fold = dsum[k0:].copy(), dsum[k0:].copy()
    #_dsum_fold[1:],_dwgt_fold[1:] = n.average(dsum_fold, axis=0), n.average(dwgt_fold, axis=0)
    #print dsum.shape, _dsum_fold.shape
    #dsum_fold,dwgt_fold = _dsum_fold,_dwgt_fold
    dsum_fold = n.concatenate([[dsum[k0]], dsum_fold])
    dwgt_fold = n.concatenate([[dwgt[k0]], dwgt_fold])
    pk_boot.append(dsum/dwgt)
    #dsum_fold = dsum[k0:].copy()
    #dwgt_fold = dwgt[k0:].copy()
    #dsum_fold[1:] = 0
    #dwgt_fold[1:] = 0
    #dsum_pos,dwgt_pos = dsum[k0+1:].copy(), dwgt[k0+1:].copy() 
    ##dsum_neg,dwgt_neg = dsum[k0-1:0:-1].copy(), dwgt[k0-1:0:-1].copy() # for even # of channels
    #dsum_neg,dwgt_neg = dsum[k0-1::-1].copy(), dwgt[k0-1::-1].copy() # for odd # of channels
    #for h in xrange(2): # bootstrap over which half of the spectrum (or both) are used
    #    h = random.randint(0,1)
    #    dsum_fold[1:] += [dsum_pos, dsum_neg][h]
    #    dwgt_fold[1:] += [dwgt_pos, dwgt_neg][h]
    pk_fold_boot.append(dsum_fold / dwgt_fold)
pk_boot = n.array(pk_boot).T
pk_fold_boot = n.array(pk_fold_boot).T

print 'Sorting bootstraps'
pk = n.average(pk_boot, axis=1)
pk_fold = n.average(pk_fold_boot, axis=1)
# this is excluding imag component in noise estimate `
pk_boot = n.sort(pk_boot.real, axis=1) # dropping imag component here
pk_fold_boot = n.sort(pk_fold_boot.real, axis=1) # dropping imag component here
if True:
    print 'Deriving errors from histogram'
    up_thresh = int(n.around(0.975 * pk_boot.shape[1])) # 2 sigma, single tail
    # important to only include real component in estimation of error
    err = (pk_boot[:,up_thresh] - pk.real) / 2 # effective "1 sigma" derived from actual 2 sigma
    up_thresh_fold = int(n.around(0.975 * pk_fold_boot.shape[1])) # 2 sigma, single tail
    err_fold = (pk_fold_boot[:,up_thresh_fold] - pk_fold.real) / 2 # effective "1 sigma" derived from actual 2 sigma
    for k_,pk_,err_ in zip(kpl, pk/1e6, 2*err/1e6):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (k_, pk_+err_,pk_,err_)
    print '-'*70
    for k_,pk_,err_ in zip(kpl[k0:], pk_fold/1e6, 2*err_fold/1e6):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (k_, pk_+err_,pk_,err_)
else:
    err = n.std(pk_boot, axis=1)
    err_fold = n.std(pk_fold_boot, axis=1)

if False:
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
    p.show()

if opts.nocov:
    print 'Writing nopspec.npz'
    n.savez('nocov_pspec.npz', kpl=kpl, pk=pk, err=err, pk_fold=pk_fold, err_fold=err_fold, cmd=cmd)
else:
    print 'Writing pspec.npz'
    n.savez(OutDIR+fname+'.npz', kpl=kpl, pk=pk, err=err, pk_fold=pk_fold, err_fold=err_fold, cmd=cmd)
    
    

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

