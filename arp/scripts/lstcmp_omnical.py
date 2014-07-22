#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a
import capo as C
import sys

def lstgridwgts(dlst, lstres):
    return n.array([-dlst/lstres, 1-n.abs(dlst/lstres), dlst/lstres]).clip(0,1).transpose()

MASKSUN = True

SEP = 'sep34'

lstres = 2 * n.pi * 42.8 / a.const.sidereal_day

if MASKSUN:
    aa = a.cal.get_aa('psa6240_v003', n.array([.15]))
    sun = a.cal.get_catalog('psa6240_v003', ['Sun'])['Sun']

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
    lst_gw = lstgridwgts(npz['lsts'] - n.array([lsts[i] for i in lst_i]), lstres)
    print len(npz['times'])
    wgt = n.where(n.isnan(npz['chi2_lin']), 0, 1./npz['chi2_lin'])
    if MASKSUN:
        for i,t in enumerate(npz['times']):
            aa.set_jultime(t)
            sun.compute(aa)
            if sun.alt > 0: wgt[i] *= 0
    gsum,gwgt = 0, 0
    info = dict(npz)
    #for sep in npz.files:
    for sep in [SEP]:
        if not sep.startswith('sep'): continue
        dat = n.where(wgt > 0, npz[sep], 0)
        lst_dat = []
        # interpolate gridded lst data to actual lst of measured data
        for i,L in enumerate(lst_i):
            d1,d2,d3 = lstnpz[sep][max(L-1,0)], lstnpz[sep][L], lstnpz[sep][min(L+1,len(lsts)-1)]
            gw1,gw2,gw3 = lst_gw[i]
            lst_dat.append(d1*gw1+d2*gw2+d3*gw3)
        lst_dat = n.array(lst_dat)
        lst_dwgt = n.array([lstnpz['wgt'][i] for i in lst_i]) # thos one for first version of lstbin
        #lst_dwgt = n.array([lstnpz['wgt_'+sep][i] for i in lst_i]) # this one for next version of lstbin
        lst_dsum = lst_dat * lst_dwgt
        
        g2 = dat * wgt * n.conj(lst_dsum) * lst_dwgt
        w2 = wgt * n.conj(lst_dsum) * lst_dsum
        gsum = n.sum(g2, axis=-1) # average across the band
        gwgt = n.sum(w2, axis=-1) # average across the band
        sep_gt = gsum/gwgt
        sep_gt.shape = (sep_gt.size, 1)
        d = lst_dsum / n.where(lst_dwgt > 0, lst_dwgt, 1)
        res = n.where(wgt > 0, dat/sep_gt - d, 0)
        sig = n.sqrt(n.median(n.abs(res)**2))
        sep_wgt = n.where(n.abs(res) > 3*sig, 0, 1)
        info['gain_'+sep] = sep_gt.flatten()
        info['wgt_'+sep] = sep_wgt
        
        #if sep == SEP:
        if False:
            gt = n.where(gwgt > 0, gsum / gwgt, 0)
            gt.shape = (gt.size, 1)
            #res1.append(n.where(wgt > 0, dat-d, 0))
            #res2.append(n.where(wgt > 0, dat/gt - d, 0))
            #p.plot(npz['lsts'], dat[:,140], '+')
            #p.plot(npz['lsts'], lst_dat[:,140], '.-')
            #p.plot(npz['lsts'], dat[:,140]/gt[:,0], '.')
            #plt_y = dat[:,140]
            #plt_x = [npz['lsts'], lsts]
            #p.plot(plt_x[0], plt_y, '+')
            #p.plot(plt_x[1], plt_y / gt[:,0], '.')
            p.subplot(121); C.arp.waterfall(dat/gt, mode='lin', mx=1000, drng=1000)
            p.colorbar(shrink=.5)
            p.subplot(122); C.arp.waterfall(res_flg, mode='lin', mx=100, drng=100)
            p.colorbar(shrink=.5)
            p.show()
    outfile = filename.split('.npz')[-2]+'_v2.npz'
    print 'Writing', outfile
    n.savez(outfile, **info)
    #C.arp.waterfall(osum/owgt, mode='lin')
    #p.colorbar()
    #p.show()
    #gt = n.where(gwgt > 0, gsum / gwgt, 0)
    #gt.shape = (gt.size, 1)
    #g_vs_t.append(gt)
    #jds.append(npz['times'])
    #chi2.append(npz['chi2_lin'])

#p.show()

#g_vs_t = n.concatenate(g_vs_t, axis=0)
#jds = n.concatenate(jds, axis=0)
#chi2 = n.concatenate(chi2, axis=0)
#res1 = n.concatenate(res1, axis=0)
#res2 = n.concatenate(res2, axis=0)
#
#p.plot(jds, n.abs(g_vs_t[:,0]))
#p.plot(jds, n.angle(g_vs_t[:,0]))
#
#p.show()
#
#p.subplot(121)
#C.arp.waterfall(res1, mode='lin', mx=200, drng=200)
#p.colorbar(shrink=.5)
#
#p.subplot(122)
#C.arp.waterfall(res2, mode='lin', mx=200, drng=200)
#p.colorbar(shrink=.5)
#
#p.show()
#
##import IPython; IPython.embed()
#import IPython.Shell; IPython.Shell.IPShellEmbed('')()


