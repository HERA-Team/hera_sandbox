#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, cal=True, src=True)
o.add_option('--bm', dest='beam', action='store_true',
    help='Divide out by model beam response.')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
srclist,cutoff,catalogs, = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)
src = cat.values()[0]
del(uv)

spec, swgt = 0, 0
for filename in args:
    uvofile = filename + 'C'
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, ants=opts.ant) 
    #uv.select('auto', -1, -1, include=False)
    data,mask,times = {}, {}, []
    print '    Reading data...'
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or times[-1] != t:
            times.append(t)
            aa.set_jultime(t)
            cat.compute(aa)
            s_eq = cat.get_crds('eq', ncrd=3)
            aa.sim_cache(s_eq)
        bl, pol = a.miriad.ij2bl(i,j), uv['pol']
        try: d = aa.phs2src(d, src, i, j)
        except(a.phs.PointingError):continue
        d /= src.get_jys()
        if opts.beam:
            d /= aa.bm_response(i,j,pol=a.miriad.pol2str[pol]).squeeze()
        if not pol in data: data[pol], mask[pol] = {}, {}
        if not bl in data[pol]: data[pol][bl], mask[pol][bl] = {}, {}
        data[pol][bl][t] = n.where(f, 0, d)
        mask[pol][bl][t] = f

    print '    Calibrating...'
    data_times, gains, wgts = None, {}, {}
    for pol in data:
      for bl in data[pol]:
        if data_times is None:
            data_times = data[pol][bl].keys()
            data_times.sort()
        data[pol][bl] = n.array([data[pol][bl][t] for t in data_times])
        mask[pol][bl] = n.array([mask[pol][bl][t] for t in data_times])
    if data_times is None or len(data_times) == 0:
        print '    No data to calibrate'
        continue
    for pol in data:
      for bl1 in data[pol]:
        i1,j1 = a.miriad.bl2ij(bl1)
        for bl2 in data[pol]:
            i2,j2 = a.miriad.bl2ij(bl2)
            if   i1 == i2 and j1 != j2:
                i,bl3 = i1,a.miriad.ij2bl(j1,j2)
                dbl1,dbl2 = data[pol][bl1], n.conj(data[pol][bl2])
                if j1 < j2: dbl3 = n.conj(data[pol][bl3])
                else: dbl3 = data[pol][bl3]
                #if i == 0: print (i1,j1), (i2,j2), (j1,j2), dbl1[10][600] * dbl2[10][600] / dbl3[10][600]
            elif i1 == j2 and j1 != i2:
                i,bl3 = i1,a.miriad.ij2bl(j1,i2)
                dbl1,dbl2 = data[pol][bl1], data[pol][bl2]
                if j1 < i2: dbl3 = n.conj(data[pol][bl3])
                else: dbl3 = data[pol][bl3]
            elif j1 == i2 and i1 != j2:
                i,bl3 = j1,a.miriad.ij2bl(i1,j2)
                dbl1,dbl2 = n.conj(data[pol][bl1]), n.conj(data[pol][bl2])
                if i1 < j2: dbl3 = n.conj(data[pol][bl3])
                else: dbl3 = data[pol][bl3]
            elif j1 == j2 and i1 != i2:
                i,bl3 = j1,a.miriad.ij2bl(i1,i2)
                dbl1,dbl2 = n.conj(data[pol][bl1]), data[pol][bl2]
                if i1 < i2: dbl3 = n.conj(data[pol][bl3])
                else: dbl3 = data[pol][bl3]
            else: continue
            if not pol in gains: gains[pol],wgts[pol] = {},{}
            if not i in gains[pol]: gains[pol][i],wgts[pol][i] = 0, 0
            f = (mask[pol][bl1] + mask[pol][bl2] + mask[pol][bl3]) > 0
            gains[pol][i] += n.where(f, 0., dbl1 * dbl2 / dbl3)
            wgts[pol][i] += n.where(f, 0., 1.)

    for pol in gains:
      for i in gains[pol]:
        #valid = wgts[pol][i] > 0
        valid = wgts[pol][i] >= wgts[pol][i].max() / 2
        gains[pol][i] = n.sqrt(n.where(valid, gains[pol][i] / wgts[pol][i], 0)).clip(5e-4,2e-2)
        wgts[pol][i] = n.where(valid, 0, 1)
    def mfunc(uv, preamble, data, flags):
        crd,t,(i,j) = preamble
        pol = uv['pol']
        try: cnt = data_times.index(t)
        except(ValueError): return preamble, None, None
        try:
            g = gains[pol][i][cnt] * n.conj(gains[pol][j][cnt])
            w = wgts[pol][i][cnt] | wgts[pol][j][cnt]
        except(KeyError): g,w = n.ones_like(data), n.zeros_like(flags)
        #if i != j:
        if False:
            data = n.where(w, 0, data / g)
            flags |= w
        else:
            data = n.where(w, 0, g * 1e6)
            flags = w
        return preamble, data, flags
        
    del(uv)
    print '    Writing output file...'
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='SELFCAL:\n')
