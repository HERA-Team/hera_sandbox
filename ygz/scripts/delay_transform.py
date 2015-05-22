__author__ = 'yunfanzhang'

import numpy as n

def nu2tau(uv1,datnu):
    datanu = datnu
    for S in n.array(range(uv1['nchan'])):
        #nu = uv1['sfreq']+S*uv1['sdf']*1.E9
        if datnu[S] == "--": datanu[S] = 0.
    datatau = n.fft.fft(datanu)
    taulist = n.fft.fftfreq(uv1['nchan'],uv1['sdf'])
    return taulist, datatau
