#! /usr/bin/env python
import aipy as a, capo as C, pylab as p, numpy as n
import sys

#ants = [1,104,106]
ants = [0,62,104,96]
#ants = [104]
antstr = ','.join(['%d_%d' % (i,i) for i in ants])
thresh = .5
tol = 1e-9
div = False
c0,c1 = 140,930
#c0,c1 = 0,1023

times,dat,flg = C.arp.get_dict_of_uv_data(sys.argv[1:], antstr=antstr, polstr='xx', verbose=True)

#colors = 'krcg'
colors = [''] * 10
g0, g1, g2 = {}, {}, {}
w0, w1, w2 = {}, {}, {}
for i,ant in enumerate(ants):
    bl = a.miriad.ij2bl(ant,ant)
    fqs = n.linspace(.1,.2,dat[bl]['xx'].shape[1])[c0:c1]
    tau = n.fft.fftfreq(fqs.size, fqs[1]-fqs[0])
    d,f = dat[bl]['xx'][:,c0:c1], flg[bl]['xx'][:,c0:c1]
    ntimes,nfreqs = d.shape
    print ant, d.shape
    w = n.logical_not(f).astype(n.float)
    d *= w
    #C.arp.waterfall(d, drng=1); p.show()
    wt = w.sum(axis=1); wt.shape = (-1,1)
    dt = n.sum(d, axis=1); dt.shape = (-1,1)
    ok = wt > thresh * wt.max()
    dt,wt = n.where(ok,dt/wt,1), n.where(ok,1,0)
    d,w = d / dt * wt, w * wt
    print n.average(w)
    #p.plot(dt); p.show()
    #C.arp.waterfall(d, drng=1); p.show()
    d,w = d.sum(axis=0), w.sum(axis=0); d /= w.max(); w /= w.max()
    ok = w > thresh
    g0[ant],w0[ant] = n.where(ok, d/w, 0), n.where(ok,1,0)
    p.subplot(121); p.plot(fqs, g0[ant], colors[i], label=str(ant), alpha=.5)
    window = a.dsp.gen_window(d.size, 'blackman-harris')
    dw,ww = d * window, w * window
    #dw,ww = g0[ant] * window, w0[ant] * window
    _dw,_ww = n.fft.ifft(dw), n.fft.ifft(ww)
    gain = a.img.beam_gain(_ww)
    mdl,info = a.deconv.clean(_dw, _ww, tol=tol, stop_if_div=div)
    mdl += info['res'] / gain
    mdl /= mdl[0]
    p.subplot(122); p.semilogy(tau[:tau.size/2], n.abs(mdl[:tau.size/2]), colors[i], label=str(ant), alpha=.5)
    sync = (fqs/.150)**-4.5
    w2[ant] = w0[ant]
    g2[ant] = n.where(w2[ant] > 0, g0[ant]/sync, 0)
    p.subplot(121); p.plot(fqs, g2[ant], colors[i], label='%d/S'%(ant))
    dw,ww = g2[ant] * window, w2[ant] * window
    _dw,_ww = n.fft.ifft(dw), n.fft.ifft(ww)
    gain = a.img.beam_gain(_ww)
    mdl,info = a.deconv.clean(_dw, _ww, tol=tol, stop_if_div=div)
    mdl += info['res'] / gain
    mdl /= mdl[0]
    p.subplot(122); p.semilogy(tau[:tau.size/2], n.abs(mdl[:tau.size/2]), colors[i], label='%d/S'%(ant))
    mdl[:16] = 0
    mdl[-15:] = 0
    p.subplot(121); p.plot(fqs, n.fft.fft(mdl), colors[i], label='%d/S'%(ant))
    p.subplot(122)
p.ylim(1e-6,1); p.xlim(0, 500); p.grid()
ns = 2./n.sqrt(100e6/1024/2*10.7*ntimes) # 2 in numerator for autovariance, 2 w/ nchan for real fft
print ns
p.plot(tau[:tau.size/2], ns*n.ones_like(tau[:tau.size/2]), 'k--')
p.xlabel('Delay (ns)')
p.ylabel('Power')
p.subplot(121)
p.xlabel('Frequency (GHz)')
p.ylabel('Power')
p.legend()
p.show()
    

