import numpy as np
from aipy import dsp,deconv
from pylab import *

c = 0.3 #m/ns
twopi = 2.*np.pi

#  _   _ _   _ _ _ _         _____       _ _
# | | | | |_(_) |_) |_ _   _|  _  \  ___| | |_
# | | | |  _| | | |  _| | | |  _  //  _ \ |  _|
# | |_| | |_| | | | |_| |_| | |_) \   __/ | |_
#  \____/\__|_|_|_|\__|__,  |_____/\____|_|\__|
#                      |____/
#

def nu2lam2(nu): return (c/nu)**2
def lam22nu(lam2): return (c/np.sqrt(lam2))
def nu_lam_swap(x): return c/x
def better_guess_l2(nu):
    dnu = np.mean(np.diff(nu))
    return (c/nu)**2 * (1. + (dnu/(4.*nu)))
def l0(nu): return np.sum(better_guess_l2(nu))/np.where(nu,1,0).sum()


#  _   _       _         ___
# | \ / | __ _| | _ __  |  _|____   ___  ___
# |  V  |/ _` | | /' _ \ \ \|  _ \ / _ \/ __|
# | |_| | (_| |   \  __/__\ | |_) |  __/ (__
# |_| |_|\__,_|_|\_\___|____/  __/ \___|\___|
#                           |_|

def gen_rm_spec(nu,rm,amp=1.):
    return amp*np.exp(2.j*rm*nu2lam2(nu))

def gen_QUspec(nu,rm,dly,amp):
    """
    input:  nu=frequencies to use
            rm=list of rotation measures
            dly=list of delays
            amp=amplitude of each 'source'
    """
    QiU = np.zeros(len(nu),dtype=np.complex)
    try:
        for i,r in enumerate(rm):
            QiU += np.exp(-1j*twopi*nu*dly[i])*gen_rm_spec(nu,r,amp[i])      
    except(TypeError): 
        QiU = np.exp(-1j*twopi*nu*dly)*gen_rm_spec(nu,rm,amp)      
    return QiU.real,QiU.imag

def dlam2(nu,dnu):
    return -2. * (c/nu)**2 * (dnu/nu)

def gen_rm_samples(nu):
    N = float(len(nu))
    j = np.arange(N)
    j -= N/2.
    L0 = l0(nu)
    nu0,dnu = np.mean(nu),nu[1]-nu[0]
    RMmax = (np.pi*nu0)/(4.*L0*dnu)
    return np.linspace(-0.5*RMmax,0.5*RMmax,N)

def Lam2Measure(nu,dnu):
    return np.abs(dlam2(nu,dnu))

#  ____        _    _
# |  _  \  ___| |__(_)_ __
# | |_) | / _ \ '_ \ | '_  \
# |  _  /|  __/ |_)| | | | |
# |_| \_\ \___|_,__/_|_| |_|
#

def rebin_nu2lam2(nu,f_nu,bin=100):
    l2 = nu2lam2(nu)
    hist1,bins = np.histogram(l2, bins=bin, weights=f_nu)
    hist2,bins = np.histogram(l2, bins=bin)
    l2 = .5 * (bins[1:] + bins[:-1])
    return l2, hist1 / np.where(hist2 == 0, 1., hist2)


#  ____  ____ _____
# |  _  \  __|_   _|
# | | | | |_   | |
# | |_| |  _|  | |
# |____ /_|    |_|
#

def RMTmat(nu,window='hamming'):
    N = len(nu)
    dnu = nu[1]-nu[0]
    W = np.indices(np.array((N,N)))
    RMs = gen_rm_samples(nu)
    L2 = better_guess_l2(nu)
    wgt = dsp.gen_window(N,window)
    W = np.exp(-2.j*RMs[W[1]]*L2[W[0]]) * Lam2Measure(nu[W[0]],dnu) * wgt[W[0]] 
    return RMs,W.T
def RMT(spec,W): return np.dot(W,spec)

#  _ ____  ____ _____
# (_)  _  \  __|_   _|
# | | | | | |_   | |
# | | |_| |  _|  | |
# |_|____ /_|    |_|
#

def iRMTmat(nu,window='hamming'):
    N = len(nu)
    dnu = nu[1]-nu[0]
    W = np.indices(np.array((N,N)))
    RMs = gen_rm_samples(nu)
    L2 = better_guess_l2(nu)
    wgt = dsp.gen_window(N,window)
    W = np.exp(2.j*RMs[W[0]]*L2[W[1]]) / wgt[W[1]]
    W *= (RMs[-1]-RMs[0])/(np.pi*N)
    return W.T
def iRMT(spec,W): return np.dot(W,spec)

#  ____  _   _  ____ _
# |  _  \ \ / |/ ___\ | ___  ___ _ _ ___
# | |_) |  V  | |   | |/ _ \/ _ ' | '_  \
# |  _  / |_| | |___| |  __/ (_|  | | | |
# |_| \_\_| |_|\____/_|\___|\___,_|_| |_|
#

def RMclean1(fq,spec,gain=0.1,tol=1e-3,stop_if_div=False,
    maxiter=10000,verbose=True,pos_def=False,window='hamming'):
    
    RMs,W = RMTmat(fq,window)
    rm_spec = RMT(spec,W)
    ker_spec = RMT(np.ones_like(spec),W)
    
    mod,info= deconv.clean(rm_spec,ker_spec,gain=gain,maxiter=maxiter,tol=tol,
        stop_if_div=stop_if_div,verbose=verbose)
    W = iRMTmat(fq,window)
    mod_spec=iRMT(mod,W)
    res_spec=iRMT(info['res'],W)

    return mod_spec,res_spec

def gen_RMdly_mat(nu,spec,window='hamming'):
    N = len(nu)
    RMs = gen_rm_samples(nu)
    wgt = dsp.gen_window(N,window)
    mat = np.zeros((N,N),dtype=np.complex)
    for i,RM in enumerate(RMs):
        mat_ift = spec*np.conjugate(gen_rm_spec(nu,RM))
        mat[:,i] = np.fft.fftshift(np.fft.ifft(mat_ift*wgt))
    return mat

def RMclean2_iter(fq,spec,RM0,DLY0,AMP,gain):
    Q,U = gen_QUspec(fq,RM0,DLY0,gain*AMP)
    return Q+1.j*U

def RMclean2(fq, QiU, gain=0.1, tol=1e-2, stop_if_div=True,window='hamming'):
    spec = QiU.copy()
    N = len(fq)
    mod = np.zeros(N,dtype=np.complex)

    dly = np.fft.fftshift(np.fft.fftfreq(N,fq[1]-fq[0]))
    rms = gen_rm_samples(fq)

    #First iteration
    W = gen_RMdly_mat(fq,spec,window=window)
    ind = np.where(np.abs(W)**2 == np.max(np.abs(W)**2))
    max0 = W[ind][0]
    max = max0
    rms_now = np.std(np.abs(spec))
    rms_last = 10.*rms_now
    rm,dl = rms[ind[1][0]],dly[ind[0][0]]
    score = 1. 
    print 0,rm,dl,score

    #the clean loop:
    i = 1
    if stop_if_div==True:
        def keep_going(max,rms_now,rms_last):
            if np.abs(max)/np.abs(max0) >= tol and rms_now < rms_last: return True
            else: return False
    else: 
        def keep_going(max,rms_now,rms_last):
            if np.abs(max)/np.abs(max0) >= tol: return True
            else: return False

    while(keep_going(max,rms_now,rms_last)): 
        mod1 = RMclean2_iter(fq,spec,rm,dl,max,gain)
        mod += mod1
        spec -= mod1
        W = gen_RMdly_mat(fq,spec)
        ind = np.where(np.abs(W) == np.max(np.abs(W)))
        score = np.abs(W[ind])/np.abs(max0)
        rm,dl = rms[ind[1]],dly[ind[0]]
        max = W[ind]
        rms_last = rms_now
        rms_now = np.std(np.abs(spec))
        print i,rm,dl,score
        i += 1
    return mod,spec

#  _   __  ____
# | | | _ \  __|
# | |_|  _/  _|
# |___|_| |_|
#

def LPF(fq,QiU,cutoff,window='hamming'):
    rms,W = RMTmat(fq,window=window)
    spec_rm = RMT(QiU,W)
    
    filter = np.where(np.abs(rms) <= cutoff,1.,0.)
    spec_rm *= filter

    W = iRMTmat(fq,window=window)
    return iRMT(spec_rm,W)

#   ___  _   _
#  / _ \| |_| |__   ___ _ __
# | | | |  _|  _ \ / _ \ '__|
# | |_| | |_| | | |  __/ |
#  \___/ \__|_| |_|\___|_|
#

def template_spectrum(nu,a):
    spec = nu + (np.pi*np.mean(nu)/l0(nu))*(c/nu)**2 
    spec = np.exp(-2.j*np.pi*a*spec)
    return spec

