#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys

# [6,11,12,13,21,24,25,26,28,30,32,35,44,47,48,49,50,52,64,84,85,99]]
#sep = 'sep34'
sep = 'sep24'

MX = 10.

d,w = [],[]
for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    d.append(npz[sep])
    try: _wgt = n.where(n.isnan(npz['chi2_lin']), 0, 1./npz['chi2_lin'])
    except(KeyError): _wgt = 1.
    wgt = npz['wgt_'+sep]
    if 'wgt_'+sep in npz.files:
        if wgt.max() == 1:
            w.append(wgt * _wgt)
        else: w.append(wgt)
    else: w.append(n.ones(d[-1].shape, dtype=n.float))
    #p.plot(n.sqrt(n.median(n.abs(npz[sep])**2, axis=0)))
d,w = n.concatenate(d),n.concatenate(w)
#w = n.where(n.isnan(d), 0, w)
#w = n.where(w < 1e-8, 0, w)
#
#d = n.where(w > 0, d, 0)
#print n.any(n.isnan(d))
#
#sig = n.sqrt(n.median(n.abs(d)**2, axis=0)); sig.shape = (1,sig.size)
#p.plot(sig[0]); p.show()
#w = n.where(n.abs(d) > 3*sig, 0, w)
#d = n.where(w > 0, d, 0)


#p.subplot(131); C.arp.waterfall(d, mx=MX, drng=MX); p.colorbar(shrink=.5)
p.subplot(131); C.arp.waterfall(d, mode='lin', mx=MX, drng=MX); p.colorbar(shrink=.5)
p.subplot(132); C.arp.waterfall(d, mode='phs'); p.colorbar(shrink=.5)
p.subplot(133); C.arp.waterfall(w, mode='log'); p.colorbar(shrink=.5)
p.show()

#import IPython.Shell; IPython.Shell.IPShellEmbed('')()

#dcpy = d.copy()
d = d[:,90:129]
w = n.where(d == 0, 0, 1.)
window = a.dsp.gen_window(d.shape[1], 'blackman-harris');
window.shape = (1,window.size)
_d = n.fft.ifft(d*window)
_w = n.fft.ifft(w*window)

p.subplot(121); C.arp.waterfall(_d)
p.subplot(122); C.arp.waterfall(_w)
p.show()

area = n.ones(d.shape[1], dtype=n.int)
#area[25:175] = 0
#area[13:190] = 0
#area[30:-31] = 0
area[9:-8] = 0

for i in range(d.shape[0]):
    print i
    g = n.sqrt(n.average(w[i]**2))
    if g == 0: continue
    _dcl,info = a.deconv.clean(_d[i], _w[i], tol=1e-9, area=area, stop_if_div=False, maxiter=100)
    dmdl = n.fft.fft(_dcl)
    d[i] -= dmdl * w[i]

'''
xtalk = n.sum(d, axis=0) / n.sum(w, axis=0)
xtalk.shape = (1,xtalk.size)
#dx = n.where(w > 0, d-xtalk, 0)
dx = n.where(w > 0, d-xtalk*w, 0)
'''

#C.arp.waterfall(n.concatenate([dcpy, d, dx]))
p.subplot(121); C.arp.waterfall(d, mx=MX, drng=MX, mode='lin')
p.colorbar(shrink=.5)
p.subplot(122); C.arp.waterfall(d, mode='phs')
p.colorbar(shrink=.5)
p.show()
