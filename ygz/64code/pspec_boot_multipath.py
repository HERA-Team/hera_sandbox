#! /usr/bin/env python
import aipy as a, numpy as np, pylab as p
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
    f = np.load(filename)
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
    #err_vs_t[path].append(np.average(f['err_vs_t'][:,120:141], axis=1))
    #temp_noise_var[path].append(np.average(f['temp_noise_var'][:,120:141], axis=1))
    nocov_vs_t[path].append(f['nocov_vs_t'])

paths = pk_vs_t.keys()
k0 = np.abs(kpl).argmin()

##############################################################################

#import IPython; IPythonp.embed()
for i,path in enumerate(paths):
    #first make all times same length
    minntimes = np.amin([arr.shape[-1] for arr in pk_vs_t[path]])
    print minntimes
    for i, arr in enumerate(pk_vs_t[path]):
        pk_vs_t[path][i] = arr[:,:minntimes]
        #import IPython; IPythonp.embed()
        nocov_vs_t[path][i] = nocov_vs_t[path][i][:,:minntimes]
        #import IPython; IPythonp.embed()
        try:
            temp_noise_dic[path][i] = temp_noise_dic[path][i][:,:minntimes]
        except:
            pass
    pk_vs_t[path] = np.array(pk_vs_t[path])
    nocov_vs_t[path] = np.array(nocov_vs_t[path])
    temp_noise_dic[path] = np.array(temp_noise_dic[path])
    #import IPython; IPythonp.embed()

pk_2d = pk_vs_t.values() # (bltype/path,bootstraps/numfilesinpath,kpls/nchan,times)
nocov_2d = nocov_vs_t.values() # (bltype,bootstraps,kpls,times), T averaged over all bls
temp_noise_var = temp_noise_dic.values()
npaths = len(pk_2d)
# try: # put everything into same numpy array if number of time samples (last axis) line up

#     pk_2d = np.stack(pk_2d) # (bltype/path,bootstraps/numfilesinpath,kpls/nchan,times)
#     nocov_2d = np.stack(nocov_2d) # (bltype,bootstraps,kpls,times), T averaged over all bls
#     temp_noise_var = np.stack(temp_noise_var)
# except:
#     pass


#import IPython; IPythonp.embed()
#average over files (bootstrap axis) in each path
avg_pk_2d = [np.average(pk_2d[i], axis=0) for i, k in enumerate(pk_2d)]
wgts = [np.ones_like(avg_pk_2d[i]) for i,av in enumerate(avg_pk_2d)]


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


#if True: # plot some stuff
if False: # plot some stuff
    plt1 = int(np.sqrt(len(paths)))
    plt2 = int(np.ceil(len(paths)/float(plt1)))

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
        #C.arp.waterfall(np.cumsum(avg_pk_2d[cnt]*wgts[cnt],axis=1)/np.cumsum(wgts[cnt],axis=1), mx=10, drng=4)
        C.arp.waterfall(np.cumsum(avg_pk_2d[cnt]*wgts[cnt],axis=1)/np.cumsum(wgts[cnt],axis=1), mode='real', mx=5e7, drng=1e8)
        p.colorbar(shrink=.5) 
    p.subplot(plt2,plt1,1); p.title('Weighted Power Spectrum [mK$^2$]'); p.show()

#print avg_pk_2d.shape, wgts.shape
for i in range(len(wgts)): print 'avg_pk_2d[i].shape, wgts[i].shape',avg_pk_2d[i].shape, wgts[i].shape


########################################################################################################
#(bltype/path,bootstraps/numfilesinpath,kpls/nchan,times) to (path,kpls/nchan,bootstraps*times)
for i in xrange(npaths):
    ntimes = pk_2d[i].shape[-1]  #now its time*bls
    nboots = pk_2d[i].shape[0]
    print 'path, ntimes, nboots =  ', i, ntimes, nboots
    #print (pk_2d[i].shape[0], nboots, ntimes), pk_2d[i].shape
    
    pk_2d[i] = pk_2d[i].transpose(1,0,2).reshape((pk_2d[i].shape[1], nboots*ntimes))

pk_2d = np.concatenate(tuple(pk_2d),axis=-1)
wgts = np.concatenate(tuple(wgts),axis=-1)
print 'pk_2d.shape',pk_2d.shape
print 'wgts.shape',wgts.shape
########################################################################################################
#ntimes = pk_2d.shape[-1] / 2

pk_boot = []
pk_fold_boot = []
for boot in xrange(NBOOT):
    if boot % 20 == 0: print boot
    #dsum,dwgt = 0, 0
    dsum,dwgt = [],[]
    dsum_fold,dwgt_fold = [], []
    for t in xrange(ntimes):
        t = np.random.randint(0, wgts.shape[-1])
        b = np.random.randint(0, pk_2d.shape[-1])
        #dsum += pk_2d[b,:,t] * wgts[:,t]
        #dwgt += wgts[:,t]
        dsum += [pk_2d[:,b] * wgts[:,t]]
        dwgt += [wgts[:,t]]
    for t in xrange(2*ntimes):
        t = np.random.randint(0, wgts.shape[-1])
        b = np.random.randint(0, pk_2d.shape[-1])
        h = np.random.randint(0,1)
        if h == 0:
            dsum_fold += [pk_2d[k0+1:,b] * wgts[k0+1:,t]]
            dwgt_fold += [wgts[k0+1:,t]]
        else:
            dsum_fold += [pk_2d[k0-1::-1,b] * wgts[k0-1::-1,t]]
            dwgt_fold += [wgts[k0-1::-1,t]]
    if MEDIAN:
        dsum,dwgt = np.median(dsum, axis=0), np.median(dwgt, axis=0)
        dsum_fold,dwgt_fold = np.median(dsum_fold, axis=0), np.median(dwgt_fold, axis=0)
    else:
        dsum,dwgt = np.average(dsum,axis=0), np.average(dwgt,axis=0)

    dsum_fold = np.concatenate([[dsum[k0]], dsum_fold])
    dwgt_fold = np.concatenate([[dwgt[k0]], dwgt_fold])
    pk_boot.append(dsum/dwgt)

    pk_fold_boot.append(dsum_fold / dwgt_fold)
pk_boot = np.array(pk_boot).T
pk_fold_boot = np.array(pk_fold_boot).T

print 'Sorting bootstraps'
pk = np.average(pk_boot, axis=1)
pk_fold = np.average(pk_fold_boot, axis=1)
# this is excluding imag component in noise estimate `
pk_boot = np.sort(pk_boot.real, axis=1) # dropping imag component here
pk_fold_boot = np.sort(pk_fold_boot.real, axis=1) # dropping imag component here
if True:
    print 'Deriving errors from histogram'
    up_thresh = int(np.around(0.975 * pk_boot.shape[1])) # 2 sigma, single tail
    # important to only include real component in estimation of error
    err = (pk_boot[:,up_thresh] - pk.real) / 2 # effective "1 sigma" derived from actual 2 sigma
    up_thresh_fold = int(np.around(0.975 * pk_fold_boot.shape[1])) # 2 sigma, single tail
    err_fold = (pk_fold_boot[:,up_thresh_fold] - pk_fold.real) / 2 # effective "1 sigma" derived from actual 2 sigma
    for k_,pk_,err_ in zip(kpl, pk/1e6, 2*err/1e6):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (k_, pk_+err_,pk_,err_)
    print '-'*70
    for k_,pk_,err_ in zip(kpl[k0:], pk_fold/1e6, 2*err_fold/1e6):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (k_, pk_+err_,pk_,err_)
else:
    err = np.std(pk_boot, axis=1)
    err_fold = np.std(pk_fold_boot, axis=1)

if False:
    #p.subplot(221); C.arp.waterfall(pk_boot, mx=9, drng=3); p.colorbar(shrink=.5)
    p.subplot(221); p.plot(pk_boot)
    p.subplot(222)
    p.plot(pk)
    p.plot(np.median(pk_boot,axis=1))
    p.plot(pk+2*err)
    #p.subplot(223); C.arp.waterfall(pk_fold_boot, mx=9, drng=3); p.colorbar(shrink=.5)
    p.subplot(223); p.plot(pk_fold_boot)
    p.subplot(224)
    p.plot(pk_fold)
    p.plot(np.median(pk_fold_boot,axis=1))
    p.plot(pk_fold+2*err_fold)
    p.show()

if opts.nocov:
    print 'Writing nopspec.npz'
    np.savez('nocov_pspec.npz', kpl=kpl, pk=pk, err=err, pk_fold=pk_fold, err_fold=err_fold, cmd=cmd)
else:
    print 'Writing pspec.npz'
    np.savez(OutDIR+fname+'.npz', kpl=kpl, pk=pk, err=err, pk_fold=pk_fold, err_fold=err_fold, cmd=cmd)


