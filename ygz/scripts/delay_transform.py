__author__ = 'yunfanzhang'

import numpy as n

def nu2tau(datnu):
    datanu = datnu
    datanu = datanu.filled(datanu.mean())
        #nu = uv1['sfreq']+S*uv1['sdf']*1.E9
    #print datanu
    datatau = n.fft.ifft(datanu)
    datatau = n.fft.fftshift(datatau)
    #print datatau

    return datatau
