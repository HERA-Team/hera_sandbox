#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import dspec, glob

rfi_file = n.load('rfi_flags.npz')
freqs = rfi_file['freqs']
#freqs = n.arange(.1,.2,.1/1024)
rfi = rfi_file['flags'].astype(n.float)
rfi = 1 - rfi/rfi.max()

print 'Reading sync_spec.npz'
sync_spec = n.load('sync_spec.npz')
umags = sync_spec['umags']
sync_spec = sync_spec['spec']
print umags

print 'Reading psrc_spec_5jy.npz'
psrc_spec = n.load('psrc_spec_5jy.npz')
psrc_spec = psrc_spec['spec']

WINDOW = 'blackman-harris'
fq0 = .150
B = .008
umag150 = 16.
Tsum,Twgt = 0, 0
#sampling = n.ones_like(freqs)
window = a.dsp.gen_window(freqs.size, window=WINDOW)
sdf = freqs[1] - freqs[0]
sampling = rfi
         
_window = n.fft.ifft(window * sampling)
bm_fqs = freqs.clip(.120,.190)

npz_dat = {'umags':umags, 'ks':[], 'fqs':[], 'spec':[]}

kwargs = {'cen_fqs':n.arange(.120,.190, B/2),'B':B, 'ntaps':3, 'window':WINDOW, 'bm_fqs':bm_fqs}
Dsum,Dwgt = {},{}
for i,umag150 in enumerate(umags):
    bl_ns = umag150 / .150
    jy_spec_i = sync_spec[i] + psrc_spec[i]
    if False:
        p.plot(freqs, sync_spec)
        p.plot(freqs, psrc_spec)
        p.plot(freqs, jy_spec_i)
        p.show()
    Trms_list = {}
    for j in xrange(1):
        d_jys, sampling = jy_spec_i.copy(), rfi.copy()
        if True: # Add noise
            noise = n.random.normal(scale=.05*(freqs/.150)**-2.5,size=freqs.size)
            noise = noise * n.exp(2j*n.pi*n.random.uniform(size=freqs.size).astype(n.complex))
            d_jys += noise
        d_jys *= sampling
        d_mdl, d_res = 0, d_jys.copy()
        for bl_frac in range(1,5):
            print bl_frac
            U,L = dspec.sky_dly_thresh(bl_ns, sdf, d_jys.size, max_bl_frac=bl_frac)
            d_res_mdl, d_res = dspec.wideband_dspec(d_res, sampling, U, L, tol=1e-9, window='blackman-harris')
            d_mdl += d_res_mdl
        if False: # Clean again w/o the window
            d_res_mdl, d_res = dspec.wideband_dspec(d_res, sampling, d_res.size/2, d_res.size/2, tol=1e-9)
            d_mdl += d_res_mdl
        if False: # Uniformize sampling by adding noise
            noise = n.random.normal(scale=.15*(freqs/.150)**-2.5,size=freqs.size)
            d_res = n.where(sampling < .5, noise, d_res/sampling)
        if False: # Plot some stuff
            p.subplot(121)
            p.plot(freqs, d_mdl,'g')
            p.plot(freqs, d_res,'r')
            p.plot(freqs, d_jys,'k')
            p.ylim(-1,1)
            p.subplot(122)
            tau = n.fft.fftfreq(freqs.size, freqs[1] - freqs[0])
            window = a.dsp.gen_window(freqs.size, window='blackman-harris')
            p.semilogy(tau, n.abs(n.fft.ifft(d_mdl*window)), 'g')
            p.semilogy(tau, n.abs(n.fft.ifft(d_res*window)), 'r')
            p.semilogy(tau, n.abs(n.fft.ifft(d_jys*window)), 'k')
            p.xlim(-1500,1500)
            p.ylim(1e-4,1e2)
            p.show()
        Tres,ks = C.pspec.Trms_vs_fq(freqs, d_res, **kwargs)
        Tmdl,ks = C.pspec.Trms_vs_fq(freqs, d_mdl, **kwargs)
        W,ks = C.pspec.Trms_vs_fq(freqs, sampling, **kwargs)
        Trms = {}
        for fq in Tres:
            _cl_Trms, info = a.deconv.clean(Tres[fq], W[fq], tol=1e-12, maxiter=100, 
                stop_if_div=False, verbose=False) 
            print int(n.around(fq*1e3)), info['term']
            Trms[fq] = Tmdl[fq] + _cl_Trms + info['res']
            if False: # Plot some stuff
                p.semilogy(ks[fq][1], n.abs(W[fq]), 'c')
                p.semilogy(ks[fq][1], n.abs(Tmdl[fq]), 'b')
                p.semilogy(ks[fq][1], n.abs(Tres[fq]), 'g')
                p.semilogy(ks[fq][1], n.abs(Trms[fq]), 'k')
                p.show()
            Trms_list[fq] = Trms_list.get(fq,[]) + [Trms[fq]]
    for fq in Trms_list:
        Di,Wi = C.pspec.k3pk_from_Trms(Trms_list[fq], [1.]*len(Trms_list[fq]), k=ks[fq][0], B=B, fq=fq)
        Dsum[fq] = Dsum.get(fq,0) + Di
        Dwgt[fq] = Dwgt.get(fq,0) + Wi

    k_img, D_img = [], []
    fq_img = []
    fqs = Dsum.keys(); fqs.sort()
    for fq in fqs:
        if True: # Use log binning
            #bins = n.arange(-3.8,3.2,1)
            bins = n.arange(-3.8,3.2,0.5)
            D = Dsum[fq] / Dwgt[fq]
            kbin,D = C.pspec.rebin_log(ks[fq][0], D, bin=bins)
            valid = n.where(n.abs(D) == 0, 0, 1)
            #kbin = kbin.compress(valid)
            D = n.where(valid, D, 1e6)
            k_img.append(kbin)
            D_img.append(D)
            fq_img.append(fq * n.ones_like(kbin))
        else:
            k_img.append(ks[fq][0])
            D_img.append(Dsum[fq] / Dwgt[fq])
            fq_img.append(fq * n.ones_like(k_img[-1]))

    # Make waterfall k3pk_vs_k_vs_fq plot
    k_img, D_img = n.array(k_img), n.array(D_img)
    fq_img = n.array(fq_img)
    #p.contourf(n.log10(k_img), fq_img, n.log10(n.abs(D_img.real).clip(1e0,1e6)), V=n.arange(0,6,.5))
    D_img = n.abs(D_img.real)
    p.contourf(n.log10(k_img), fq_img, n.log10(D_img.clip(1e1,1e6)), levels=n.arange(1,6.5,.5))
    npz_dat['fqs'].append(fq_img)
    npz_dat['ks'].append(k_img)
    npz_dat['spec'].append(D_img)
    #p.colorbar(shrink=.5)
    p.xlim(-1.5, 0.31)
    p.xlabel(r'${\rm log}_{10}\,k\,[h\,{\rm Mpc}^{-1}]$', size=14)
    p.ylabel('Frequency (GHz)', size=14)
    p.show()

print 'Writing fg_vs_umag_vs_fq.npz'
npz_dat['fqs'] = n.array(npz_dat['fqs'])
npz_dat['ks'] = n.array(npz_dat['ks'])
npz_dat['spec'] = n.array(npz_dat['spec'])
n.savez('fg_vs_umag_vs_fq.npz', **npz_dat)
