#! /usr/bin/env python
import numpy as n, aipy as a
import sys

seps = [6,11,12,13,21,24,25,26,28,30,32,35,44,47,48,49,50,52,64,84,85,99]
seps = ['sep%d'%i for i in seps]

def lstgridwgts(dlst, lstres):
    return n.array([-dlst/lstres, 1-n.abs(dlst/lstres), dlst/lstres]).clip(0,1).transpose()

MASKSUN = False # already done in lstcmp_omnical
INTERPOLATE = True

lstres = 2 * n.pi * 42.8 / a.const.sidereal_day
lstsum,lstwgt = {}, {}

if MASKSUN:
    aa = a.cal.get_aa('psa6240_v003', n.array([.15]))
    sun = a.cal.get_catalog('psa6240_v003', ['Sun'])['Sun']

for filename in sys.argv[1:]:
    print 'Reading', filename
    npz = n.load(filename)
    lsts = n.around(npz['lsts'] / lstres) * lstres
    wgt = n.where(n.isnan(npz['chi2_lin']), 0, 1./npz['chi2_lin'])
    lst_gw = lstgridwgts(lsts-npz['lsts'], lstres)
    if MASKSUN:
        for i,t in enumerate(npz['times']):
            aa.set_jultime(t)
            sun.compute(aa)
            if sun.alt > 0: wgt[i] *= 0
    for sep in seps:
        if not sep.startswith('sep'): continue
        if not lstsum.has_key(sep):
            lstsum[sep],lstwgt[sep] = {}, {}
        if 'gain_'+sep in npz.files:
            #gt = npz['gain_'+sep]
            gt = npz['gain_'+sep].clip(.5,1.5)
            gt.shape = (gt.size,1)
        else: gt = 1
        if 'wgt_'+sep in npz.files: sep_wgt = wgt * npz['wgt_'+sep]
        else: sep_wgt = wgt
        dat = n.where(sep_wgt > 0, npz[sep]/gt, 0)
        sep_wgt = n.where(n.isnan(dat), 0, sep_wgt)
        dat = n.where(sep_wgt > 0, dat, 0)
        print sep, n.any(n.isnan(dat))
        if INTERPOLATE:
            _dat,_wgt = [], []
            for i in range(len(lsts)):
                d1,d2,d3 = dat[max(i-1,0)], dat[i], dat[min(i+1,len(lsts)-1)]
                w1,w2,w3 = sep_wgt[max(i-1,0)], sep_wgt[i], sep_wgt[min(i+1,len(lsts)-1)]
                gw1,gw2,gw3 = lst_gw[i]
                _dat.append(d1*w1*gw1+d2*w2*gw2+d3*w3*gw3)
                _wgt.append(w1*gw1+w2*gw2+w3*gw3)
            _wgt= n.array(_wgt)
            dat = n.array(_dat) / n.where(_wgt > 0, _wgt, 1)
        else: _wgt = sep_wgt

        # done in wideband_dly_omnical
        #if True:
        #    window = a.dsp.gen_window(dat[:,14:186].shape[1], 'blackman-harris'); 
        #    window.shape = (1,window.size)
        #    w = n.where(_wgt > 0, 1., 0)
        #    _d = n.fft.ifft(dat[:,14:186]*window)
        #    _w = n.fft.ifft(w[:,14:186]*window)

        #    area = n.ones(dat[:,14:186].shape[1], dtype=n.int)
        #    area[19:-18] = 0 # XXX gotta make this bl dependent

        #    for i in range(dat.shape[0]):
        #        g = n.sqrt(n.average(w[i]**2))
        #        if g < .5 :
        #            dat[i], w[i], _wgt[i] = 0, 0, 0
        #            continue
        #        _dcl,info = a.deconv.clean(_d[i], _w[i], tol=1e-9, area=area, stop_if_div=False, maxiter=100)
        #        dmdl = n.fft.fft(_dcl)
        #        dat[i,14:186] -= dmdl * w[i,14:186]

        #    xtalk = n.sum(dat, axis=0) / n.sum(w, axis=0)
        #    xtalk.shape = (1,xtalk.size)
        #    dat = n.where(w > 0, dat-xtalk, 0)

        
        for i,lst in enumerate(lsts):
            lstsum[sep][lst] = lstsum[sep].get(lst,0) + dat[i] * _wgt[i]
            lstwgt[sep][lst] = lstwgt[sep].get(lst,0) + _wgt[i]

lsts = lstwgt.values()[0].keys(); lsts.sort()
data = {}

for sep in lstsum.keys():
    print sep
    d,w = [],[]
    for lst in lsts:
        d.append(lstsum[sep].pop(lst))
        w.append(lstwgt[sep].pop(lst))
    del(lstsum[sep]); del(lstwgt[sep])
    d,w = n.array(d), n.array(w)
    data[sep] = d / n.where(w > 0, w, 1)
    data['wgt_'+sep] = w

filename = 'omnical_lstbin_v2_yy.npz'
print 'Writing', filename
n.savez(filename, lsts=n.array(lsts), **data)
