__author__  =  'yunfanzhang'
import aipy as a, numpy as n
from scipy import interpolate

#computes the beam squared
def beam_real(ant, ntop, shape0=None, pol='x', sq=True):
        Nfreq = len(ant.bm_response(ntop,pol=pol))
        bmp_list = []
        for nf in range(Nfreq):
            bm1x = ant.bm_response(ntop,pol=pol)[nf]
            bm2x = ant.bm_response(ntop,pol=pol)[nf]
            bm = bm1x*n.conj(bm2x)
            bmsq = bm*bm
            bmp = bm
            if sq: bmp = bmsq
            if shape0!=None: bmp.shape = shape0
            bmp_list.append(bmp)
        return bmp_list

def OmP(ant, ntop, pol, sq=False):
    bm = ant.bm_response(ntop,pol=pol)
    bm = bm*n.conj(bm)
    dl = abs(ntop[0][2]-ntop[0][1])
    if sq: bm = bm*bm
    #print bm.shape
    #print dl
    return n.sum(bm*dl*dl)

#computes the fourier transform of give beam pattern bmp
def beam_fourier(bmp, dreal, nreal, dl=0.005):
        fbm = n.fft.fft2(bmp)
        frequv = n.fft.fftfreq(nreal,d=dreal)
        freqk = frequv*2*n.pi
        #fbmamp = fbm.real*dl*dl
        fbmamp = n.abs(fbm)*dl*dl
        #fbmamp = n.abs(fbm)
        freq = frequv
        fbmamp = n.array(n.fft.fftshift(fbmamp))
        freq = n.array(n.fft.fftshift(freq))
        return freq, fbmamp

def beam_interpol(freq, fbmamp,kind='cubic'):
    return interpolate.interp2d(freq, freq, fbmamp, kind=kind)

#Interpolates for the overlap of two baselines given a (u,v) coordinate,i.e. the Fourier transform of the beam squared
def get_overlap(f, u, v, Diag=False):
        if len(f(u,v)) > 1:
            if Diag: return n.diagonal(f(u,v))
            else: return f(u,v)
        else:
            return f(u,v)[0]