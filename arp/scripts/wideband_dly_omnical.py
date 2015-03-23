#! /usr/bin/env python
import numpy as n, aipy as a
import sys

for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    info = dict(npz)
    _wgt = n.where(n.isnan(npz['chi2_lin']), 0, 1./npz['chi2_lin'])
    #for sep in npz.files:
    #for sep in ['sep34']:
    for sep in ['sep%d'%i for i in [6,11,12,13,21,24,25,26,28,30,32,35,44,47,48,49,50,52,64,84,85,99]]:
        if not sep.startswith('sep'): continue
        wgt = npz['wgt_'+sep]
        wgt = n.where(n.isnan(wgt), 0, wgt)
        wgt = n.where(n.isnan(npz[sep]), 0, wgt)
        dat = n.where(wgt > 0, npz[sep], 0)
        print sep, n.any(n.isnan(dat))

        # XXX edges of band hardcoded here
        window = a.dsp.gen_window(dat[:,14:186].shape[1], 'blackman-harris'); 
        window.shape = (1,window.size)
        w = n.where(wgt > 0, 1., 0)
        _d = n.fft.ifft(dat[:,14:186]*window)
        _w = n.fft.ifft(w[:,14:186]*window)

        area = n.ones(dat[:,14:186].shape[1], dtype=n.int)
        area[19:-18] = 0 # XXX gotta make this bl dependent

        for i in range(dat.shape[0]):
            g = n.sqrt(n.average(w[i]**2))
            if g < .5 :
                dat[i], wgt[i] = 0, 0
                continue
            _dcl,junk = a.deconv.clean(_d[i], _w[i], tol=1e-9, area=area, stop_if_div=False, maxiter=100)
            dmdl = n.fft.fft(_dcl)
            dat[i,14:186] -= dmdl * w[i,14:186]

        #xtalk = n.sum(dat, axis=0) / n.sum(w, axis=0)
        xtalk = n.sum(dat*_wgt * wgt, axis=0) / n.sum(_wgt * wgt, axis=0)
        xtalk.shape = (1,xtalk.size)
        dat = n.where(w > 0, dat-xtalk, 0)
        info[sep] = dat
        info['wgt_'+sep] = wgt
    outfile = filename.split('.npz')[-2]+'B.npz'
    print 'Writing', outfile
    n.savez(outfile, **info)

