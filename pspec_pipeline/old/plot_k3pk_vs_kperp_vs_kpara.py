#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C

freqs = n.linspace(.1,.2,1024)
WINDOW = 'blackman-harris'
fq0 = .150
B = .008
umag150 = 16.
bl_ns = umag150 / .150

#def sim_synchrotron(fq_GHz, amp_Jy_150=400e3, index=-2.5, umag_index=-3., umag_150=40., umag_150_0=3.):
#    spec = amp_Jy_150 * (u150/u150_0 * freqs/mfreq)**index_vs_umag * (freqs / mfreq)**index

def sim_src_Jy(fq_GHz, amp_Jy_150=.1, index=-1., dly_ns=100., phs_offset=1.):
    spec = amp_Jy_150 * (fq_GHz / .150)**index
    phs = phs_offset * n.exp(2j*n.pi * fq_GHz.astype(n.complex64) * dly_ns)
    return spec * phs
    
def sim_srcs(fq_GHz, cnt_100Jy_per_sr_Jy=.1, cnt_index=-2., avg_index=-1., std_index=.25, 
        lo_cutoff_Jy=.3, hi_cutoff_Jy=100., bl_ns=100.):
    bins = 10**n.arange(n.log10(hi_cutoff_Jy), n.log10(lo_cutoff_Jy), -.3)
    spec = n.zeros_like(fq_GHz, dtype=n.complex64)
    for i, flx in enumerate(bins[:-1]):
        flx_interval = bins[i] - bins[i+1]
        flx += flx_interval/2
        cnt = int(n.around(cnt_100Jy_per_sr_Jy * flx_interval * 2*n.pi * (flx / 100.)**cnt_index))
        # I'm ignoring beam weighting on sky and geometric weight toward low delays by just selecting
        # random delays and random phase offsets
        for j in xrange(cnt):
            index = avg_index + n.random.normal(scale=std_index)
            dly = n.random.uniform(-bl_ns,bl_ns)
            phs = n.exp(2j*n.pi*n.random.uniform())
            spec += sim_src_Jy(fq_GHz, flx, index=index, dly_ns=dly, phs_offset=phs)
    return spec
   
Tsum,Twgt = 0, 0
sampling = n.ones_like(freqs)
#sampling[500:502] = 0
#sampling[520:525] = 0
#sampling[530:531] = 0
window = a.dsp.gen_window(freqs.size, window=WINDOW) * sampling
_window = n.fft.ifft(window)
bm_fqs = freqs.clip(.120,.190)
#etas = C.pspec.f2eta(freqs)
#kwargs = {'cen_fqs':[fq0], 'B':B, 'ntaps':3, 'window':WINDOW, 'bm_fqs':bm_fqs}
kwargs = {'cen_fqs':n.arange(.130,.180,.005),'B':B, 'ntaps':3, 'window':WINDOW, 'bm_fqs':bm_fqs}
Dsum,Dwgt = {},{}

# Make k3pk vs kpara vs kperp plot
kwargs = {'cen_fqs':[fq0], 'B':B, 'ntaps':3, 'window':WINDOW, 'bm_fqs':bm_fqs}
Dsum, Dwgt = {}, {}
k_vs_umag = {}
umags = n.arange(20, 300, 40)
sampling = n.ones_like(freqs)
for umag in umags:
    bl_ns = umag / .150
    Trms_list = []
    jy_spec_i = sim_srcs(freqs, bl_ns=bl_ns, lo_cutoff_Jy=.1, hi_cutoff_Jy=5., std_index=.25)
    for i in xrange(5):
        noise = n.random.normal(scale=.05*(freqs/.150)**-2.5,size=freqs.size)
        jy_spec = jy_spec_i + noise
        jy_spec *= sampling
        _jy_spec = n.fft.ifft(jy_spec*window)
        _cl_jy_spec, info = a.deconv.clean(_jy_spec, _window, tol=1e-9, stop_if_div=False)
        jy_spec2 = jy_spec - n.fft.fft(_cl_jy_spec) * sampling
        jy_spec_mdl = n.fft.fft(_cl_jy_spec)

        Trms1,ks = C.pspec.Trms_vs_fq(freqs, jy_spec, umag150=umag, **kwargs)
        Trms2,ks = C.pspec.Trms_vs_fq(freqs, jy_spec2, umag150=umag, **kwargs)
        Tmdl ,ks = C.pspec.Trms_vs_fq(freqs, jy_spec_mdl, umag150=umag, **kwargs)
        W,ks = C.pspec.Trms_vs_fq(freqs, sampling, umag150=umag, **kwargs)
        Trms5 = {}
        _cl_Trms, info = a.deconv.clean(Trms2[fq0], W[fq0], tol=1e-9, maxiter=100, 
            stop_if_div=False, verbose=False) 
        print i, int(umag), info['term']
        Trms5[fq0] = Tmdl[fq0] + _cl_Trms + info['res']
        Trms_list.append(Trms5[fq0])
    Dsum[umag] = C.pspec.k3pk_from_Trms(Trms_list, k=ks[fq0][0], B=B, fq=fq0)
    Dwgt[umag] = 1
    k_vs_umag[umag] = ks[fq0]

kpara, kperp, D = [], [], []
for umag in umags:
    k,_kpara,_kperp = k_vs_umag[umag]
    if False:
        kbin,Dbin = C.pspec.rebin_log(k, Dsum[umag] / Dwgt[umag], nbins=20)
        #valid = n.where(n.abs(Dbin) == 0, 0, 1)
        #kbin = kbin.compress(valid)
        Dbin = n.where(n.abs(Dbin) == 0, 1e6, Dbin)
        _kpara = n.sqrt(kbin**2 - _kperp**2)
        kpara.append(_kpara)
        kperp.append(_kperp * n.ones_like(_kpara))
        D.append(Dbin)
    else:
        kpara.append(_kpara)
        kperp.append(_kperp * n.ones_like(_kpara))
        D.append(Dsum[umag] / Dwgt[umag])
kpara,kperp,D = n.abs(kpara).clip(1e-3,n.Inf), n.abs(kperp).clip(1e-3,n.Inf), n.array(D)
print kpara, kperp, D

#p.contourf(n.log10(kperp), n.log10(kpara), n.log10(n.abs(D.real).clip(1e0,1e6)), 20)
#p.contourf(n.log10(kperp), n.log10(kpara), n.log10(n.abs(D.real).clip(1e0,1e6)), levels=n.arange(0,6.5,.5))
#p.contourf(kperp, n.log10(kpara), n.log10(n.abs(D.real).clip(1e0,1e6)), levels=n.arange(0,6.5,.5))
p.contourf(kperp, kpara, n.log10(n.abs(D.real).clip(1e0,1e6)), levels=n.arange(0,6.5,.5))
#p.colorbar(shrink=.5)
p.ylim(-1, 0.31)
p.xlabel(r'$k_\perp\,[h\,{\rm Mpc}^{-1}]$', size=14)
p.ylabel(r'$k_\parallel\,[h\,{\rm Mpc}^{-1}]$', size=14)
p.show()
