'''Units in mK, GHz, Mpc unless stated otherwise'''
import numpy as n, aipy as a

F21 = 1.42040575177 # GHz
LST_RES = 2*n.pi/24
UV_RES = 1.5
# Cosmological conversions
def f2z(fq):
    '''Convert frequency (GHz) to redshift for 21cm line.'''
    return (F21 / fq - 1)
def z2f(z):
    '''Convert redshift to frequency (GHz) for 21cm line.'''
    return F21 / (z+1)
def dL_df(z):
    '''[h^-1 Mpc]/GHz, from Furlanetto et al. (2006)'''
    return (1.7 / 0.1) * ((1+z) / 10.)**.5 * 1e3
def dL_dth(z):
    '''[h^-1 Mpc]/radian, from Furlanetto et al. (2006)'''
    return 1.9 * (1./a.const.arcmin) * ((1+z) / 10.)**.2
def dk_deta(z):
    '''2pi * [h Mpc^-1] / [GHz^-1]'''
    return 2*n.pi / dL_df(z) 
def dk_du(z):
    '''2pi * [h Mpc^-1] / [wavelengths]'''
    return 2*n.pi / (dL_dth(z) * (n.pi / 2))
def X2Y(z):
    '''[h^-3 Mpc^3] / [str * GHz] scalar conversion between observing and cosmological coordinates'''
    return dL_dth(z)**2 * dL_df(z)

# Reionization models
def k3pk_21cm(k):
    '''Return peak \Delta_{21}^2(k), in mK**2, approximation to Furlanetto et al. (2006)'''
    return n.where(k < 0.1, 1e4 * k**2, 1e2)
def dk3pk_21cm(k, dkx, dky, dkz):
    '''Return peak 21cm pspec, integrated over dk**3 box size'''
    return k3pk_21cm(k) / (4*n.pi * k**3) * dkx * dky * dkz

# Temperature conversion
DEFAULT_BEAM_POLY = [ -1.55740671e+09,  1.14162351e+09, -2.80887022e+08,  9.86929340e+06, 7.80672834e+06, -1.55085596e+06,  1.20087809e+05, -3.47520109e+03]

def jy2T(f, bm_poly=DEFAULT_BEAM_POLY):
    '''Return [mK] / [Jy] for a beam size vs. frequency (in GHz) defined by the
    polynomial bm_poly.  Default is a poly fit to the PAPER primary beam.'''
    lam = a.const.c / (f * 1e9)
    bm = n.polyval(bm_poly, f)
    return 1e-23 * lam**2 / (2 * a.const.k * bm) * 1e3
def k3pk_from_Trms(Trms, k=.3, fq=.150, B=.001, bm_poly=DEFAULT_BEAM_POLY):
    z = f2z(fq)
    bm = n.polyval(bm_poly, fq)
    return X2Y(z) * bm * B * k**3 / (2*n.pi)**2 * Trms**2
def k3pk_sense_vs_t(t, k=.3, fq=.150, B=.001, bm_poly=DEFAULT_BEAM_POLY, Tsys=500e3):
    Trms = Tsys / n.sqrt(2*(B*1e9)*t)
    return k3pk_from_Trms(Trms, k=k, fq=fq, B=B, bm_poly=bm_poly)

# Misc helper functions
def uv2bin(u,v,lst,uv_res=UV_RES, lst_res=LST_RES):
    return (int(n.around(u / uv_res) + 4096) * 8192 + int(n.around(v / uv_res) + 4096)) * 8192 + lst/lst_res
def bin2uv(bin, uv_res=UV_RES, lst_res=LST_RES):
    v = ((bin/8192) % 8192 - 4096) * float(uv_res)
    u = (bin / 8192**2 - 4096) * float(uv_res)
    lst = (bin % 8192) * float(lst_res)
    return u,v, lst
def rebin_log(x, y, nbins=10):
    '''For y=f(x), bin x into log_10 bins, and average y over
    these bin sizes.'''
    logx = n.log10(n.abs(x))
    hist1,bins = n.histogram(logx, bins=nbins, weights=y)
    hist2,bins = n.histogram(logx, bins=nbins)
    logx = .5 * (bins[1:] + bins[:-1])
    return 10**logx, hist1 / n.where(hist2 == 0, 1., hist2)
def f2eta(f):
    '''Convert an array of frequencies to an array of etas (freq^-1) 
    corresponding to the bins that an FFT generates.'''
    return n.fft.fftfreq(f.shape[-1], f[1]-f[0])
