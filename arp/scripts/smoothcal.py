#! /usr/bin/env python

import numpy as np, aipy, capo, pylab as plt, sys, glob

aa = aipy.cal.get_aa('psa6622_v003', np.array([.15]))

SEPS = [
    (0,103) , #  1
    (1,4) ,   #  2
    (0,101) , #  3
    (0,62) ,  #  4
    (0,100) , #  5
    (1,13) ,  #  6
    (1,70) ,  #  7
    (1,56) ,  #  8
    (1,71) ,  #  9
    (1,59) ,  # 10
    (0,97) ,  # 11
    (12,43) , # 12
    (9,71) ,  # 13
    (9,59) ,  # 14
    (57,64) , # 15
]

CONJ = [
    (0,103) , #  1
    (1,4) ,   #  2
    (0,101) , #  3
    (0,62) ,  #  4
    (0,100) , #  5
    (0,97) ,  # 11
    (12,43) , # 12
    (57,64) ] # 15


POL = 'xx'
files = glob.glob('zen.2456715.*.%s.npz' % POL)
meta, gains, data, xtalk = capo.omni.from_npz(files[0], verbose=True)
    
NITER = 30
fq = np.linspace(.1,.2,203)
window = aipy.dsp.gen_window(fq.size, 'hamming')
tau = np.fft.fftfreq(fq.size, fq[1]-fq[0])
tau_sh = np.fft.fftshift(tau)
bp = np.ones(fq.size, dtype=np.complex)
for cnt in xrange(NITER):
    bp_sum, bp_wgt = 0., 0.
    #_bp = np.fft.ifft(bp, axis=-1)
    #for (i,j) in data[POL]:
    for (i,j) in SEPS:
        bl_vec = aa[j] - aa[i] # bl_vec in units of ns
        tau_mx = np.sqrt(np.dot(bl_vec,bl_vec))
        #if tau_mx > 200: continue
        print cnt, (i,j), bl_vec, tau_mx
        d = data[POL][(i,j)][0].copy()
        ok = np.where(np.abs(d) > 0, 1, 0) * window
        _bp = np.fft.ifft(bp * ok, axis=-1)
        if bl_vec[1] < 0: d = d.conj()
        #d /= np.where(np.abs(bp) > 0, bp, 1)
        _d = np.fft.ifft(d * window, axis=-1)
        area = np.where(np.abs(tau) < tau_mx, 1, 0)
        mdl, info = aipy.deconv.clean(_d, _bp, area=area, tol=1e-6, stop_if_div=False, maxiter=100)
        print info['term']
        if True:#and cnt > 0:
            plt.subplot(121)
            mdl0,info0 = aipy.deconv.clean(_d, np.fft.ifft(ok), area=area, tol=1e-6, stop_if_div=False, maxiter=100)
            #plt.semilogy(tau_sh, np.fft.fftshift(np.abs(_d)))
            plt.semilogy(tau_sh, np.fft.fftshift(np.abs(mdl0)), 'r')
            plt.semilogy(tau_sh, np.fft.fftshift(np.abs(info0['res'])), 'r')
            plt.semilogy(tau_sh, np.fft.fftshift(np.abs(mdl)), 'k')
            plt.semilogy(tau_sh, np.fft.fftshift(np.abs(info['res'])), 'k')
            #plt.show()
        d_s = np.fft.fft(mdl) * bp
        d_f = np.fft.fft(info['res']) / window
        bpi = (d_f / d_s + 1)
        if True:#and cnt > 0:
            plt.subplot(122)
            #plt.plot(bp, 'r')
            plt.plot(bpi)
            #plt.show()
        # XXX think about averaging in log space
        bp_sum += bpi / tau_mx
        bp_wgt += 1. / tau_mx
        #plt.subplot(121)
        #plt.plot(bpi)
        #plt.subplot(122)
        #plt.semilogy([tau_mx,tau_mx],[1e-6,1], 'r:')
        #plt.semilogy([-tau_mx,-tau_mx],[1e-6,1], 'r:')
        #plt.semilogy(tau_sh, np.fft.fftshift(np.abs(_d)**2), 'k.')
    plt.show()
    bpi = (bp_sum / bp_wgt)**.2
    #plt.plot(bp, 'k.')
    #plt.plot(bpi, 'r.')
    bp *= bpi
    #plt.plot(bp, 'b.')

plt.plot(np.abs(bp), 'k')
plt.plot(bp.real, 'b')
plt.plot(bp.imag, 'r')
plt.show()

for (i,j) in data[POL]:
    bl_vec = aa[j] - aa[i] # bl_vec in units of ns
    tau_mx = np.sqrt(np.dot(bl_vec,bl_vec))
    print (i,j), bl_vec, tau_mx
    d = data[POL][(i,j)][0].copy()
    if bl_vec[1] < 0: d = d.conj()
    bp0 = np.where(np.abs(d) > 0, 1, 0)
    _d = np.fft.ifft(d * window, axis=-1)
    _bp0 = np.fft.ifft(bp0 * window, axis=-1)
    _bp = np.fft.ifft(bp * bp0 * window, axis=-1)
    mdl,info = aipy.deconv.clean(_d,_bp,tol=1e-9, stop_if_div=False)
    _dc = mdl + info['res']
    mdl,info = aipy.deconv.clean(_d,_bp0,tol=1e-9, stop_if_div=False)
    _dc0 = mdl + info['res']
    plt.semilogy([tau_mx,tau_mx],[1e-6,1], 'g:')
    plt.semilogy([-tau_mx,-tau_mx],[1e-6,1], 'g:')
    plt.semilogy(tau_sh, np.fft.fftshift(np.abs(_d)), 'r.')
    plt.semilogy(tau_sh, np.fft.fftshift(np.abs(_dc0)), 'b.')
    plt.semilogy(tau_sh, np.fft.fftshift(np.abs(_dc)), 'k.')
    plt.show()
    
import IPython; IPython.embed()

