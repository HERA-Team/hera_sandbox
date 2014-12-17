#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a
import capo as C
import sys

def lstgridwgts(dlst, lstres):
    return n.array([-dlst/lstres, 1-n.abs(dlst/lstres), dlst/lstres]).clip(0,1).transpose()

MASKSUN = True

SEP = 'sep6'

lstres = 2 * n.pi * 42.8 / a.const.sidereal_day

#if MASKSUN:
#    aa = a.cal.get_aa('psa6240_v003', n.array([.15]))
#    sun = a.cal.get_catalog('psa6240_v003', ['Sun'])['Sun']

lstnpz = n.load(sys.argv[1])
lsts = lstnpz['lsts']

g_vs_t = []
jds = []
chi2 = []
res1,res2 = [], []

for filename in sys.argv[2:]:
    print 'Reading', filename,
    npz = n.load(filename)
    lst_i = [n.argmin(n.abs(lsts - lst)) for lst in npz['lsts']]
    lst_dwgt = n.array([lstnpz['wgt'][i] for i in lst_i])
    lst_gw = lstgridwgts(n.array([lsts[i] for i in lst_i])-npz['lsts'], lstres)
    print len(npz['times'])
    wgt = n.where(n.isnan(npz['chi2_lin']), 0, 1./npz['chi2_lin'])
#    if MASKSUN:
#        for i,t in enumerate(npz['times']):
#            aa.set_jultime(t)
#            sun.compute(aa)
#            if sun.alt > 0: wgt[i] *= 0
    wgt[npz['sunflags']] *= 0 
    gsum,gwgt = 0, 0
    #for sep in npz.files:
    #for sep in ['sep%d'%i for i in range(10)]:
    for sep in [SEP]:
        if not sep.startswith('sep'): continue
        dat = n.where(wgt > 0, npz[sep], 0)
        _dat,_wgt = [],[] # XXX could move wgt out
        # XXX should instead interpolate lst data to measured lst
        for i in range(len(lst_i)):
            d1,d2,d3 = dat[max(i-1,0)], dat[i], dat[min(i+1,len(lst_i)-1)]
            w1,w2,w3 = wgt[max(i-1,0)], wgt[i], wgt[min(i+1,len(lst_i)-1)]
            gw1,gw2,gw3 = lst_gw[i]
            _dat.append(d1*w1*gw1+d2*w2*gw2+d3*w3*gw3)
            _wgt.append(w1*gw1+w2*gw2+w3*gw3)
        _dat = n.array(_dat)
        _wgt= n.array(_wgt)
        dat = _dat / n.where(_wgt > 0, _wgt, 1)
        #_dat,_wgt = dat,wgt
        
        lst_dsum = n.array([lstnpz[sep][i] for i in lst_i]) * lst_dwgt
        g2 = dat * _wgt * n.conj(lst_dsum) * lst_dwgt
        w2 = _wgt * n.conj(lst_dsum) * lst_dsum
        #gsum += g2
        #gwgt += w2
        gsum += n.sum(g2, axis=-1)
        gwgt += n.sum(w2, axis=-1)
        if sep == SEP:
            gt = n.where(gwgt > 0, gsum / gwgt, 0)
            gt.shape = (gt.size, 1)
            d = lst_dsum / n.where(lst_dwgt > 0, lst_dwgt, 1)
            res1.append(n.where(_wgt > 0, dat-d, 0))
            res2.append(n.where(_wgt > 0, dat/gt - d, 0))
            #plt_y = dat[:,140]
            #plt_x = [npz['lsts'], lsts]
            #p.plot(plt_x[0], plt_y, '+')
            #p.plot(plt_x[1], plt_y / gt[:,0], '.')
    gt = n.where(gwgt > 0, gsum / gwgt, 0)
    gt.shape = (gt.size, 1)
    g_vs_t.append(gt)
    jds.append(npz['times'])
    chi2.append(npz['chi2_lin'])

#p.show()

g_vs_t = n.concatenate(g_vs_t, axis=0)
jds = n.concatenate(jds, axis=0)
chi2 = n.concatenate(chi2, axis=0)
res1 = n.concatenate(res1, axis=0)
res2 = n.concatenate(res2, axis=0)

p.plot(jds, n.abs(g_vs_t[:,0]))
p.plot(jds, n.angle(g_vs_t[:,0]))

p.show()

p.subplot(121)
C.arp.waterfall(res1, mode='lin', mx=200, drng=200)
p.colorbar(shrink=.5)

p.subplot(122)
C.arp.waterfall(res2, mode='lin', mx=200, drng=200)
p.colorbar(shrink=.5)

p.show()

#import IPython; IPython.embed()
import IPython.Shell; IPython.Shell.IPShellEmbed('')()

