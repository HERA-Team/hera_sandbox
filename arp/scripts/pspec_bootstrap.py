#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

LOG = True

def noise(d):
    return n.sqrt(n.average(n.abs(n.concatenate([d.imag[:d.size/4], d.imag[3*d.size/4:]]))**2))
    #return n.sqrt(n.average(n.abs(n.concatenate([d.real[:d.size/4], d.real[3*d.size/4:]]))**2))

def plotit(kpl, d, log=LOG, colors='kb'):
    noise_est = noise(d)
    p.subplot(121)
    if log:
        #p.semilogy(kpl, d.real, 'k')
        #p.semilogy(kpl, n.abs(d.real), 'k')
        p.semilogy(kpl, d.real + 2*noise_est, colors[0])
        #p.semilogy(kpl, n.where(n.abs(n.angle(d)) < .5, d.real, 0), colors[0]+'.')
        p.semilogy(kpl, n.where(d.real > 0, d.real, 0), colors[0]+'.')
        p.semilogy(kpl, noise_est*n.ones_like(d.imag), colors[1]+':')
        p.semilogy(kpl, 2*noise_est*n.ones_like(d.imag), colors[1]+'--')
    else:
        p.plot(kpl, d.real, colors[0])
        p.plot(kpl, noise_est*n.ones_like(d.imag), colors[1]+':')
        p.plot(kpl, 2*noise_est*n.ones_like(d.imag), colors[1]+'--')
    p.subplot(122)
    k3 = n.abs(n.sqrt(kpl**2 + .01**2)**3 / (2*n.pi**2)) # XXX hardcode for |u|=16 baselines
    #k3 = n.abs(kpl**3 / (2*n.pi**2))
    p.semilogy(kpl, k3 * (d.real + 2*noise_est), colors[0])
    p.semilogy(kpl, k3 * n.abs(d.real), colors[0]+'.')
    #p.semilogy(kpl, k3 * (d.real - 2*noise_est).clip(1e-6,n.Inf), colors[0])
    p.semilogy(kpl, k3 * noise_est*n.ones_like(d.imag), colors[1]+':')
    p.semilogy(kpl, k3 * 2*noise_est*n.ones_like(d.imag), colors[1]+'--')

NBOOT = 1000
#for filename in args:
if True:
    filename = args[0]
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
    fq0 = n.average(freqs)
    etas = n.fft.fftshift(C.pspec.f2eta(freqs))
    kpl = etas * C.pspec.dk_deta(C.pspec.f2z(fq0))
    print kpl
    del(uv)
    #times, dat, flg = C.arp.get_dict_of_uv_data([filename], opts.ant, opts.pol, verbose=True)
    times, dat, flg = C.arp.get_dict_of_uv_data(args, opts.ant, opts.pol, verbose=True)
    print dat.keys()
    print len(dat.keys())
    boot = []
    for i in xrange(NBOOT):
        dsum, dwgt = 0., 0.
        print '%d/%d' % (i+1,NBOOT)
        N = len(dat.keys())
        for j in n.random.randint(0,N,N):
            k = dat.keys()[j]
            dsum += dat[k]['yy']
            dwgt += n.logical_not(flg[k]['yy']).astype(n.int)
        boot.append(dsum/dwgt)
    boot = n.array(boot)
    print boot.shape
    avg_2d = n.average(boot.real, axis=0)
    std_2d = n.std(boot.real, axis=0)
    p.subplot(131); C.arp.waterfall(avg_2d, mode='real'); p.colorbar(shrink=.5)
    p.subplot(132); C.arp.waterfall(std_2d, mode='real'); p.colorbar(shrink=.5)
    boot = []
    for i in xrange(NBOOT):
        dsum, dwgt = 0., 0.
        print '%d/%d' % (i+1,NBOOT)
        N = avg_2d.shape[0]
        for j in n.random.randint(0,N,N):
            w = 1./std_2d[j]**2
            dsum += avg_2d[j] * w
            dwgt += w
        boot.append(dsum/dwgt)
    boot = n.array(boot)
    avg_1d = n.average(boot, axis=0)
    std_1d = n.std(boot, axis=0)
    p.subplot(133)
    p.plot(kpl, avg_1d, 'k.')
    p.plot(kpl, avg_1d+2*std_1d, 'k')
    p.show()
sys.exit(0)
    
            
    #a.scripting.uv_selector(uv, opts.ant, opts.pol)
    #for (uvw,t,(i,j)),d,f in uv.all(raw=True):
    #    noise_est = noise(d)
    #    if noise_est == 0: continue
    #    print (i,j), t, noise_est / 1e6#, d.real[24:27] / 1e6
    #    dsum += d / noise_est**2
    #    dwgt += 1./ noise_est**2
    #    #plotit(kpl, d, colors='cm')
d = dsum / dwgt
nos = noise(d)
for k,pk in zip(kpl, d.real/1e6): print k, pk, '+/-', nos / 1e6
for rng in [(-0.5,-0.3),(-0.3,-0.2),(-0.2,-0.1),(0.1,0.2),(0.2,0.3),(0.3,0.5)]:
    xerr = (rng[1] - rng[0]) / 2
    rng = n.where(n.logical_and(kpl >= rng[0], kpl < rng[1]), 1, 0)
    rng_kpl = n.sum(rng * kpl) / n.sum(rng)
    rng_pk = n.sum(rng * d.real) / n.sum(rng)
    rng_nos = nos / n.sqrt(n.sum(rng))
    rng_k3 = n.abs(rng_kpl)**3 / (2*n.pi**2)
    print 'Best:', rng_kpl, rng_pk/1e6, '+/-', rng_nos/1e6
    print 'Best:', rng_kpl, rng_k3*rng_pk, '+/-', rng_k3 * rng_nos
    p.subplot(121)
    p.errorbar(rng_kpl, rng_pk+2*rng_nos, xerr=xerr, yerr=([min(4*rng_nos,rng_pk+2*rng_nos-1e-6)], [0.]), fmt='r')
    p.subplot(122)
    p.errorbar(rng_kpl, rng_k3*(rng_pk+2*rng_nos), 
        xerr=xerr, yerr=([min(rng_k3*4*rng_nos,rng_k3*(rng_pk+2*rng_nos)-1e-6)],[0.]), fmt='r')
plotit(kpl, d)
p.subplot(121)
p.gca().set_yscale('log')
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$P(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$')
p.ylim(1e5,1e9)
p.subplot(122)
p.gca().set_yscale('log')
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$k^3/2\pi\ P(k)\ [{\rm mK}^2]$')
p.ylim(1e0,1e6)
p.show()
