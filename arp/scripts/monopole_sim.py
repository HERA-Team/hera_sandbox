#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p

fqs = n.linspace(.1,.2,203)
#TSKY,TGND = 300,100 # K
TSKY,TGND = 1.,0 # K

BM_FQ = .15
#aa = a.cal.get_aa('psa898_v003', n.array([BM_FQ]))
aa = a.cal.get_aa('psa6240_v003', n.array([BM_FQ]))
#h = a.healpix.HealpixMap(nside=128)
h = a.healpix.HealpixMap(nside=64)
tx,ty,tz = h.px2crd(n.arange(h.map.size), ncrd=3)
#blx, bly, blz = aa.get_baseline(0,3,'z')
#blx, bly, blz = aa.get_baseline(0,16,'z')
blx, bly, blz = aa.get_baseline(1,4,'z')
#blx, bly, blz = aa.get_baseline(0,8,'z')
print blx, bly, blz
bl_prj = tx*blx + ty*bly + tz*blz
bm = aa[0].bm_response((tx,ty,tz),pol='x')[0]**2
bm = n.where(tz > 0, bm, 0); bm /= bm.sum()
#C.arp.plot_hmap_ortho(h, mode='real'); p.show()
spec = n.array([n.sum(TSKY*bm*n.exp(2j*n.pi*bl_prj*fq)) for fq in fqs])
window = a.dsp.gen_window(fqs.size, window='blackman-harris')
_spec = n.fft.ifft(window*spec) / n.fft.ifft(window)[0]
tau = n.fft.fftfreq(fqs.size, fqs[1]-fqs[0])
tau,_spec = n.fft.fftshift(tau), n.fft.fftshift(_spec)
p.subplot(121); p.plot(fqs, spec.real, fqs, spec.imag)
p.subplot(122); p.semilogy(tau, n.abs(_spec))
p.show()
