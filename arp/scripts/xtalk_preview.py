#! /usr/bin/env python
import aipy as a, capo as C, pylab as p, numpy as n
import sys

CH0,CH1 = 16,182
fqs = n.linspace(.1,.2,203)
jy2T = C.pspec.jy2T(fqs)
fqs = fqs[CH0:CH1]
aa = a.cal.get_aa('psa6240_v003', n.array([.15]))
t,dat,flg = C.arp.get_dict_of_uv_data(sys.argv[1:], 'cross', 'I')
window = a.dsp.gen_window(fqs.size,'blackman-harris')
norm = n.fft.ifft(window)[0]
tau = n.fft.fftfreq(fqs.size, fqs[1]-fqs[0])
tau = n.fft.fftshift(tau)

#for filename in sys.argv[1:]:
#  t,dat,flg = C.arp.get_dict_of_uv_data([filename], 'cross', 'I')
for bl in dat:
    i,j = a.miriad.bl2ij(bl)
    print i,j
    for pol in dat[bl]:
        d = n.sum(dat[bl][pol], axis=0) * jy2T
        if (aa[i] - aa[j])[1] < 0: d = d.conj()
        w = n.sum(n.logical_not(flg[bl][pol]).astype(n.int), axis=0)
        d = n.where(w[CH0:CH1] > 0, d[CH0:CH1]/w[CH0:CH1], 0)
        w = n.where(w > 0, 1, 0)
        p.subplot(121); p.plot(fqs, d)
        _d = n.fft.ifft(window*d) / norm
        _d = n.fft.fftshift(_d)
        p.subplot(122); p.plot(tau, n.abs(_d))
p.xlabel('Delay [ns]')
p.ylabel('Power [mK]')
p.show()
