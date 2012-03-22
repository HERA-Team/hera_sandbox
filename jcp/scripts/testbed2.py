#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse, os, pylab as pl

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True,cmap=True)
o.add_option('--clean', dest='clean', type='float',
    help='Deconvolve delay-domain data by the "beam response" that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts,args = o.parse_args(sys.argv[1:])

data = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, data['sdf'], data['sfreq'], data['nchan'])
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
del(data)

l = 0
m = 1

def dt(uv, p, d):
    #delay transform data
    mask = d.mask
    dint = n.ma.absolute(d.filled(0))
    dint = n.ma.masked_less_equal(dint, 0)
    din.append(n.ma.log10(dint))
    dflags = n.logical_not(d.mask).astype(n.float)
    dgain = n.sqrt(n.average(dflags**2))
    dker = n.fft.ifft(dflags)
    d = d.filled(0)
    d = n.fft.ifft(d)
    if not opts.clean is None and not n.all(d == 0):
        d, info = a.deconv.clean(d, dker, tol=opts.clean)
        d += info['res'] / dgain
    d = n.ma.array(d)
    d = n.ma.concatenate([d[d.shape[0]/2:],
                                  d[:d.shape[0]/2]], axis=0)
    d.mask = n.zeros_like(mask)
    doutt = n.ma.absolute(d.filled(0))
    doutt = n.ma.masked_less_equal(doutt, 0)
    dtout.append(n.ma.log10(doutt))
    return p, d
    #return uvt, dtout

def idt(uv, p, d):
    #inverse delay transform data
    mask = d.mask
    #dflags = n.logical_not(d.mask).astype(n.float)
    #dgain = n.sqrt(n.average(dflags**2))
    #dker = n.fft.ifft(dflags)
    #d = d.filled(0)
    d = n.fft.fft(d)
    #if not opts.clean is None and not n.all(d == 0):
    #    d, info = a.deconv.clean(d, dker, tol=opts.clean)
    #    d += info['res'] / dgain
    d = n.ma.array(d)
    #d = n.ma.concatenate([d[d.shape[0]/2:],
    #                              d[:d.shape[0]/2]], axis=0)
    d.mask = n.zeros_like(mask)
    idoutt = n.ma.absolute(d.filled(0))
    idoutt = n.ma.masked_less_equal(idoutt, 0)
    idtout.append(n.ma.log10(idoutt))
    return p, d

def phs(uv, p, d):
    #phase data source
    (uvw, t, (i, j)) = p
    #dinpt = n.ma.absolute(d.filled(0))
    #dinpt = n.ma.masked_less_equal(dinpt, 0)
    #dinp.append(n.ma.log10(dinpt))
    aa.set_jultime(t)
    cat.compute(aa)
    phses = aa.gen_phs(cat[opts.src], i, j)
    phses = n.fft.ifft(phses)
    phses = n.ma.array(phses)
    phses = n.ma.concatenate([phses[phses.shape[0]/2:],
                                  phses[:phses.shape[0]/2]], axis=0)
    #dold = d
    d = d * phses
    #print d[450:460] - dold[450:460]
    doutp = n.ma.absolute(d.filled(0))
    doutp = n.ma.masked_less_equal(doutp, 0)
    dpout.append(n.ma.log10(doutp))
    return p, d

for filename in args:
    print 'Reading', filename
    uvi = a.miriad.UV(filename)
    uvi.select('antennae', l, m)
    din = []
    dtout = []
    dpout = []
    dinp = []
    idtout = []

    print 'delay transforming...'
    uvt = a.miriad.UV('temp.uv', status='new')
    uvt.init_from_uv(uvi)
    uvt.pipe(uvi, mfunc=dt)
    del(uvt)
    uvt = a.miriad.UV('temp.uv')
    uvt.rewind()
    print 'done'
    print 'phasing data...'
    uvp = a.miriad.UV('temp2.uv', status='new')
    uvp.init_from_uv(uvt)
    uvp.pipe(uvt, mfunc=idt)
    del(uvp)
    print 'done'

    din = n.array(din)
    dtout = n.array(dtout)
    dpout = n.array(dpout)
    idtout = n.array(idtout)
    
    cmap = pl.get_cmap(opts.cmap)
    pl.subplot(311)
    pl.imshow(din, aspect='auto', cmap=cmap)
    pl.subplot(312)
    pl.imshow(idtout, aspect='auto', cmap=cmap)
    pl.subplot(313)
    pl.colorbar()
    pl.imshow(din-idtout, aspect='auto', cmap=cmap)
    pl.show()



    #uvo.init_from_uv(uvi, override=override)
    #uvo.pipe(uvi, mfunc=mfunc)
del(uvi, uvt)
os.system('rm -r temp.uv')
os.system('rm -r temp2.uv')
