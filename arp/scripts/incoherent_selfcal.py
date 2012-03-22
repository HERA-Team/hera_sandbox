#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    uvofile = filename + 'K'
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, ants=opts.ant) 
    data,wgts = {},{}
    print '    Reading data...'
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        d = n.where(f, 0, n.abs(d)**2)
        v = n.where(f, 0, 1)
        pol = uv['pol']
        if not pol in data: data[pol], wgts[pol] = {}, {}
        for ant in [i,j]:
            if not ant in data[pol]: data[pol][ant], wgts[pol][ant] = {}, {}
            data[pol][ant][t] = data[pol][ant].get(t,0) + d
            wgts[pol][ant][t] = wgts[pol][ant].get(t,0) + v
    print '    Calibrating...'
    tsum,twgt = {}, {}
    for pol in data:
        for ant in data[pol]:
            for t in data[pol][ant]:
                tsum[t] = tsum.get(t,0) + data[pol][ant][t]
                twgt[t] = twgt.get(t,0) + wgts[pol][ant][t]

    #print '    Calibrating...'
    #data_times, gains, wgts = None, {}, {}
    #for pol in data:
    #  for bl in data[pol]:
    #    if data_times is None:
    #        data_times = data[pol][bl].keys()
    #        data_times.sort()
    #    data[pol][bl] = n.array([data[pol][bl][t] for t in data_times])
    #    mask[pol][bl] = n.array([mask[pol][bl][t] for t in data_times])
    #if data_times is None or len(data_times) == 0:
    #    print '    No data to calibrate'
    #    continue
    #for pol in data:
    #  for bl1 in data[pol]:
    #    i1,j1 = a.miriad.bl2ij(bl1)
    #    for bl2 in data[pol]:
    #        i2,j2 = a.miriad.bl2ij(bl2)
    #        if   i1 == i2 and j1 != j2:
    #            i,bl3 = i1,a.miriad.ij2bl(j1,j2)
    #            dbl1,dbl2 = data[pol][bl1], n.conj(data[pol][bl2])
    #            if j1 < j2: dbl3 = n.conj(data[pol][bl3])
    #            else: dbl3 = data[pol][bl3]
    #            #if i == 0: print (i1,j1), (i2,j2), (j1,j2), dbl1[10][600] * dbl2[10][600] / dbl3[10][600]
    #        elif i1 == j2 and j1 != i2:
    #            i,bl3 = i1,a.miriad.ij2bl(j1,i2)
    #            dbl1,dbl2 = data[pol][bl1], data[pol][bl2]
    #            if j1 < i2: dbl3 = n.conj(data[pol][bl3])
    #            else: dbl3 = data[pol][bl3]
    #        elif j1 == i2 and i1 != j2:
    #            i,bl3 = j1,a.miriad.ij2bl(i1,j2)
    #            dbl1,dbl2 = n.conj(data[pol][bl1]), n.conj(data[pol][bl2])
    #            if i1 < j2: dbl3 = n.conj(data[pol][bl3])
    #            else: dbl3 = data[pol][bl3]
    #        elif j1 == j2 and i1 != i2:
    #            i,bl3 = j1,a.miriad.ij2bl(i1,i2)
    #            dbl1,dbl2 = n.conj(data[pol][bl1]), data[pol][bl2]
    #            if i1 < i2: dbl3 = n.conj(data[pol][bl3])
    #            else: dbl3 = data[pol][bl3]
    #        else: continue
    #        if not pol in gains: gains[pol],wgts[pol] = {},{}
    #        if not i in gains[pol]: gains[pol][i],wgts[pol][i] = 0, 0
    #        f = (mask[pol][bl1] + mask[pol][bl2] + mask[pol][bl3]) > 0
    #        gains[pol][i] += n.where(f, 0., dbl1 * dbl2 / dbl3)
    #        wgts[pol][i] += n.where(f, 0., 1.)

    #for pol in gains:
    #  for i in gains[pol]:
    #    #valid = wgts[pol][i] > 0
    #    valid = wgts[pol][i] >= wgts[pol][i].max() / 2
    #    gains[pol][i] = n.sqrt(n.where(valid, gains[pol][i] / wgts[pol][i], 0)).clip(5e-4,2e-2)
    #    wgts[pol][i] = n.where(valid, 0, 1)
    def mfunc(uv, preamble, dat, flg):
        crd,t,(i,j) = preamble
        pol = uv['pol']
        if i != j: return preamble, None, None
        d,w = data[pol][i][t], wgts[pol][i][t]
        f = n.where(w < w.max()/2, 1, 0)
        d = n.where(f, 0, n.sqrt(d/w))
        td,tw = tsum[t],twgt[t]
        tf = n.where(tw < tw.max()/2, 1, 0)
        td = n.where(tf, 0, n.sqrt(td/tw))
        d = n.where(tf, 0, d/td)
        return preamble, d, f
        
    del(uv)
    print '    Writing output file...'
    uvi = a.miriad.UV(filename)
    uvi.select('auto',-1,-1,include=True)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')
