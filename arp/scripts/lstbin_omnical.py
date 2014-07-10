#! /usr/bin/env python
import numpy as n, aipy as a
import sys

def lstgridwgts(dlst, lstres):
    return n.array([-dlst/lstres, 1-n.abs(dlst/lstres), dlst/lstres]).clip(0,1).transpose()

MASKSUN = True
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
    if INTERPOLATE:
        lst_gw = lstgridwgts(lsts-npz['lsts'], lstres)
        _wgt = []
        for i in range(len(lsts)):
            w1,w2,w3 = wgt[max(i-1,0)], wgt[i], wgt[min(i+1,len(lsts)-1)]
            gw1,gw2,gw3 = lst_gw[i]
            _wgt.append(w1*gw1+w2*gw2+w3*gw3)
        _wgt= n.array(_wgt)
    else: _wgt = wgt
    if MASKSUN:
        for i,t in enumerate(npz['times']):
            aa.set_jultime(t)
            sun.compute(aa)
            if sun.alt > 0: wgt[i] *= 0
    for sep in npz.files:
        if not sep.startswith('sep'): continue
        if not lstsum.has_key(sep):
            lstsum[sep],lstwgt = {}, {}
        dat = n.where(wgt > 0, npz[sep], 0)
        if INTERPOLATE:
            _dat = []
            for i in range(len(lsts)):
                d1,d2,d3 = dat[max(i-1,0)], dat[i], dat[min(i+1,len(lsts)-1)]
                w1,w2,w3 = wgt[max(i-1,0)], wgt[i], wgt[min(i+1,len(lsts)-1)]
                gw1,gw2,gw3 = lst_gw[i]
                _dat.append(d1*w1*gw1+d2*w2*gw2+d3*w3*gw3)
            dat = n.array(_dat) / n.where(_wgt > 0, _wgt, 1)
        
        for i,lst in enumerate(lsts):
            lstsum[sep][lst] = lstsum[sep].get(lst,0) + dat[i] * _wgt[i]
            lstwgt[lst] = lstwgt.get(lst,0) + _wgt[i]

lsts = lstwgt.keys(); lsts.sort()
data = {'wgt': n.array([lstwgt[lst] for lst in lsts])}

for sep in lstsum.keys():
    print sep
    d = []
    for lst in lsts:
        _d,_w = lstsum[sep].pop(lst),lstwgt[lst]
        d.append(_d / n.where(_w > 0, _w, 1))
    del(lstsum[sep])
    data[sep] = n.array(d)

filename = 'omnical_lstbin.npz'
print 'Writing', filename
n.savez(filename, lsts=n.array(lsts), **data)
