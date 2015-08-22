__author__ = 'yunfanzhang'
import aipy as a, numpy as n
from scipy import interpolate

def get_overlap(ant, ntop, shape0, pol, u, v,dreal,nreal):
        bm1x=ant.bm_response(ntop,pol='x')[0]
        bm2x=ant.bm_response(ntop,pol='x')[0]

        bm=bm1x*n.conj(bm2x)
        bmsq=(bm1x*n.conj(bm2x))*(bm1x*n.conj(bm2x))

        #Tranform the square of the beams
        bmp=bmsq
        bmp.shape=shape0

        fbm=n.fft.fft2(bmp)
        frequv=n.fft.fftfreq(nreal,d=dreal)
        freqk=frequv*2*n.pi
        fbmamp=n.log(fbm.real)
        #fbmamp=n.abs(fbm)

        freq=frequv
        fbmamp=n.fft.fftshift(fbmamp)
        freq=n.fft.fftshift(freq)

        f = interpolate.interp2d(freq, freq, fbmamp, kind='cubic')

        print f(u,v)
        return f(u,v)