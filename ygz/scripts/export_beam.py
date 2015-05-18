__author__  =  'yunfanzhang'
import aipy as a, numpy as n
from scipy import interpolate

#computes the beam squared
def beam_real(ant, ntop, shape0, pol, sq=True):
        Nfreq = len(ant.bm_response(ntop,pol=pol))
        bmp_list = []
        for nf in range(Nfreq):
            bm1x = ant.bm_response(ntop,pol=pol)[nf]
            bm2x = ant.bm_response(ntop,pol=pol)[nf]
            bm = bm1x*n.conj(bm2x)
            bmsq = bm*bm
            bmp = bm
            if sq: bmp = bmsq
            bmp.shape = shape0
            bmp_list.append(bmp)
        return bmp_list

#computes the fourier transform of give beam pattern bmp
def beam_fourier(bmp, dreal, nreal):
        fbm = n.fft.fft2(bmp)
        frequv = n.fft.fftfreq(nreal,d=dreal)
        freqk = frequv*2*n.pi
        fbmamp = fbm.real
        #fbmamp = n.abs(fbm)
        freq = frequv
        fbmamp = n.fft.fftshift(fbmamp)
        freq = n.fft.fftshift(freq)
        return freq, fbmamp

#Interpolates for the overlap of two baselines given a (u,v) coordinate
def get_overlap(freq,fbmamp, u, v, Diag=False):
        f  =  interpolate.interp2d(freq, freq, fbmamp, kind='cubic')
        if len(f(u,v)) > 1:
            if Diag: return n.diagonal(f(u,v))
            else: return f(u,v)
        else:
            return f(u,v)[0]