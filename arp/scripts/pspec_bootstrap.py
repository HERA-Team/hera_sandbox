#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

NBOOT1 = 10
BOOTLEN1 = 100
#BOOTLEN1 = 300
NBOOT2 = 10
#BOOTLEN2 = 1000
BOOTLEN2 = 100
#for filename in args:
if True:
    filename = args[0]
    uv = a.miriad.UV(filename)
    freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
    fq0 = n.average(freqs)
    etas = n.fft.fftshift(C.pspec.f2eta(freqs))
    kpl = etas * C.pspec.dk_deta(C.pspec.f2z(fq0))
    #print kpl
    del(uv)
    #times, dat, flg = C.arp.get_dict_of_uv_data([filename], opts.ant, opts.pol, verbose=True)
    times, dat, flg = C.arp.get_dict_of_uv_data(args, opts.ant, opts.pol, verbose=True)
    print dat.keys()
    N = len(dat.keys())
    for k in dat.keys():
        i,j = a.miriad.bl2ij(k)
        print (i/32,i%32), (j/32,j%32), dat[k][opts.pol].shape
    # bootstrap across baselines, derive SNR for each integration
    if N <= 1:
        avg_2d = dat.values()[0].values()[0]
        std_2d = n.ones_like(avg_2d)
    else:
      avg_2ds,std_2ds = [], []
      for h in xrange(NBOOT1):
        boot = []
        boot_dat = []
        print 'Bootstrapping over %d baseline products' % N
        for j in n.arange(0,N):
            k = dat.keys()[j]
            boot_dat.append(dat[k][opts.pol])
        for i in xrange(BOOTLEN1):
          dsum, dwgt = 0., 0.
          if i % 10 == 0: print '%d/%d-%d/%d' % (h,NBOOT1,i,BOOTLEN1)
          N = len(dat.keys())
          for j in n.random.randint(0,N,N):
              k = dat.keys()[j]
              dsum += dat[k][opts.pol]
              dwgt += n.logical_not(flg[k][opts.pol]).astype(n.int)
          boot.append(dsum/dwgt)
        boot = n.array(boot)
        print boot.shape
        boot_dat = n.array(boot_dat)
        #p.clf(); #p.plot(boot_dat[:,:,26].real.flatten(), boot_dat[:,:,26].imag.flatten(), ',')
        #p.plot(boot[:,:,26].real, boot[:,:,26].imag,','); p.show()
        avg_2d = n.average(boot, axis=0)
        #med_2d = n.median(boot.real, axis=0)
        std_2d = n.std(boot, axis=0)
        std_2d = n.where(std_2d == 0, 1e-6, 0)
        avg_2ds.append(avg_2d)
        std_2ds.append(std_2d)
      avg_2ds = n.array(avg_2ds); avg_2d = n.average(avg_2ds, axis=0)
      std_2ds = n.array(std_2ds); std_2d = n.average(std_2ds, axis=0)
    if False: # smooth out std
        std_2d_sm = n.median(std_2d, axis=-1)
        std_2d_sm.shape = std_2d_sm.shape + (1,)
        std_2d = n.ones_like(std_2d) * std_2d_sm
    p.subplot(141); C.arp.waterfall(avg_2d, mode='log', mx=9, drng=3); p.colorbar(shrink=.5)
    p.subplot(142); C.arp.waterfall(std_2d, mode='log', mx=9, drng=3); p.colorbar(shrink=.5)
    #p.subplot(143); C.arp.waterfall(avg_2d/std_2d**2, mode='log', mx=-6, drng=4); p.colorbar(shrink=.5)
    p.subplot(143); C.arp.waterfall(avg_2d/std_2d**2, mode='real', mx=1e-7, drng=2e-7); p.colorbar(shrink=.5)
    # Loop for bootstrapping weights for each integration
    avg_1ds = []
    std_1ds = []
    for h in xrange(NBOOT2):
      boot = []
      for i in xrange(BOOTLEN2):
        dsum, dwgt = 0., 0.
        if i % 10 == 0: print '%d/%d-%d/%d' % (h,NBOOT2,i,BOOTLEN2)
        N = avg_2d.shape[0]
        #p.clf(); p.plot(avg_2d[:,26].real, avg_2d[:,26].imag); p.show()
        for j in n.random.randint(0,N,N):
            #w = 1./std_2d[j]**2
            w = 1.
            dsum += avg_2d[j] * w
            dwgt += w
        boot.append(dsum/dwgt)
      boot = n.array(boot)
      #p.clf(); p.plot(avg_2d[:,26].real, avg_2d[:,26].imag, '.')
      #p.plot(boot[:,26].real, boot[:,26].imag, 'x'); p.show()
      print boot.shape
      avg_1d = n.average(boot, axis=0)
      #med_1d = n.median(boot, axis=0)
      std_1d = n.std(boot, axis=0)
      avg_1ds.append(avg_1d)
      std_1ds.append(std_1d)
    avg_1ds = n.array(avg_1ds); avg_1d = n.average(avg_1ds, axis=0)
    std_1ds = n.array(std_1ds); std_1d = n.average(std_1ds, axis=0)
    #print avg_1ds[:,26]
    #print std_1ds[:,26]
    p.subplot(144)
    p.plot(kpl, avg_1d, 'k.')
    #p.plot(kpl, med_1d, 'k+')
    p.plot(kpl, avg_1d+2*std_1d, 'k')
    filename = 'pspec.npz'
    print 'Saving', filename
    n.savez(filename, kpl=kpl, pk=avg_1d, err=std_1d)
    #for _kpl,_pk,_ns in zip(kpl, avg_1d, std_1d):
    #    print '%7.3f: (%f,%f),' % (_kpl, _pk/1e6, _ns/1e6)
    #p.show()
