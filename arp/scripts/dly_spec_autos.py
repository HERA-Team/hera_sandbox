#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys

REBIN = False

uv = a.miriad.UV(sys.argv[-1])
sdf = uv['sdf']
freqs = uv['sfreq'] + n.arange(uv['nchan'])*sdf
del(uv)
CH0,CH1 = 300,-250
#CH0,CH1 = -500,-250
freqs = freqs[CH0:CH1]
fq = n.average(freqs)
z = C.pspec.f2z(fq)

plot_data = []
for cnt, filename in enumerate(sys.argv[1:]):
  print filename
  for pcnt, pol in enumerate(['xx','yy']):
    #t, dat, flg = C.arp.get_dict_of_uv_data([filename], 'auto,-40,-55', pol)
    t, dat, flg = C.arp.get_dict_of_uv_data([filename], 'auto,-20,-27', pol)
    for acnt, bl in enumerate(dat):
      #i,i = a.miriad.bl2ij(bl)
      acnt += pcnt * 30
      print acnt
      d = dat[bl]
      d = n.average(d, axis=0)
      d_df = d.copy(); d_df[1:-1] -= (d_df[:-2] + d_df[2:])/2
      sig = n.sqrt(n.median(n.abs(d_df)**2))
      #print sig
      d = n.where(n.abs(d_df) > 10*sig, 0, d)
      d = d[CH0:CH1]
      d_df = d_df[CH0:CH1]
      w = n.where(d == 0, 0, 1.)
      if True:
          poly = n.polyfit(freqs.compress(w), n.log10(d.compress(w)), deg=6)
          d = n.where(d == 0, 0, d - 10**n.polyval(poly,freqs))
          gain = 10**n.polyval(poly,freqs) / (300e3*(freqs/.150)**-2.5 + 100)
          d /= gain # Fake putting in K units
      window = a.dsp.gen_window(d.size, window='kaiser3')
      #window = a.dsp.gen_window(d.size, window='hanning')
      #window = 1
      _d = n.fft.ifft(d * window)
      _w = n.fft.ifft(w * window)
      eta = n.fft.fftfreq(d.size, sdf)
      g = n.sqrt(n.average((w*window)**2))
      _d_c,info = a.deconv.clean(_d, _w, tol=1e-4)
      _d_c += info['res'] / g
      eta = eta[:_d_c.size/2]
      k = C.pspec.dk_deta(z) * eta
      _d_c = _d_c[:_d_c.size/2]
      B = freqs[-1] - freqs[0]
      if False:
          plot_data.append(_d_c)
      else:
          if cnt == 0: plot_data.append(_d_c)
          else: plot_data[acnt] = (plot_data[acnt] - _d_c) / n.sqrt(2)

k = n.sqrt(k**2 + .006**2)
scalar = C.pspec.k3pk_from_Trms(1, k=k, fq=fq, B=B)
plot_data = n.array(plot_data)
sum, wgt = 0, 0
for i in range(plot_data.shape[0]):
    for j in range(i+1, plot_data.shape[0]):
        sum += plot_data[i] * n.conj(plot_data[j])
        wgt += 1

p.subplot(122)
_k, _dat = k, scalar*n.abs(plot_data[0])**2
if REBIN: _k, _dat = C.pspec.rebin_log(_k, _dat, nbins=30)
p.loglog(_k, _dat/32)

p.subplot(121)
plot_data = n.log10(n.abs(plot_data)**2)
p.imshow(plot_data, vmax=3, vmin=3-3, aspect='auto',
    interpolation='nearest', origin='lower')

p.subplot(122)
_k, _dat = k, scalar*sum/wgt
if REBIN: _k, _dat = C.pspec.rebin_log(_k, _dat, nbins=30)
p.loglog(_k, n.abs(_dat))
p.show()

    
