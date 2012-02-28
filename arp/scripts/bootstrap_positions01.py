#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, chan=True, ant=True, src=True, pol=True, dec=True)
o.add_option('--clean', dest='clean', type='float', default=1e-3)
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))
uv = a.miriad.UV(args[0])
chan = a.scripting.parse_chans(opts.chan, uv['nchan'])
delays = n.arange(-.5/uv['sdf'], .5/uv['sdf'], 1/(uv['sdf']*uv['nchan']))
delays = n.concatenate([delays[delays.size/2:], delays[:delays.size/2]])
delays *= -1
delays += 55
#delays /= a.const.len_ns # cm/m
#delays /= 100

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)

data = {}

for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if n.all(f): continue
        if not data.has_key(t):
            data[t] = {'srcs':{}, 'bls':{}}
            aa.set_jultime(t)
            cat.compute(aa)
            for src in cat:
                if cat[src].alt <= 0: continue
                data[t]['srcs'][src] = cat[src].get_crds('top')
        bl = a.miriad.ij2bl(i,j)
        valid = n.logical_not(f).astype(n.float)
        gain = n.sqrt(n.average(valid**2))
        d *= valid
        ker = n.fft.ifft(valid)
        _d = n.fft.ifft(d)
        _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
        _d += info['res'] / gain
        _d = n.abs(_d)
        data[t]['bls'][bl] = _d

SIZE = 2000
RES = 10
DIM = SIZE/RES

def gen_line(slope, offset):
    dx, dy = offset
    x = n.arange(-SIZE/2, SIZE/2, RES/4)
    y = slope * x
    return x+dx, y+dy

field = {}
for t in data:
    for bl in data[t]['bls']:
        if not field.has_key(bl): field[bl] = n.zeros((DIM,DIM), dtype=n.float)
        dat = data[t]['bls'][bl]
        for src in data[t]['srcs']:
            sx,sy,sz = data[t]['srcs'][src]
            for d,dly in zip(dat,delays):
                #lin_x, lin_y = gen_line(-sx/sy, (0,dly/sy))
                crds = n.array(gen_line(-sx/sy, (0,dly/sy)))
                crds = crds.transpose()
                valid = n.where((crds**2).sum(axis=1) < (SIZE/2 - RES)**2, 1, 0)
                crds = crds.compress(valid, axis=0)
                crds = n.around(crds / RES).astype(n.int) + DIM/2
                w = n.ones(crds.shape[0], dtype=n.float)
                d = d * w
                v = field[bl][crds[:,0], crds[:,1]]
                v = n.where(v == 0, d, v)
                field[bl][crds[:,0], crds[:,1]] = n.where(v < d, v, d)
                #a.utils.add2array(field,crds,d)
                #a.utils.add2array(wgts,crds,w)

x,y = n.indices((DIM,DIM))
x,y = x-DIM/2, y-DIM/2
window = n.where(x**2 + y**2 < (DIM/2-2)**2, 1, 0)
err_width = 10/RES
err_kernel = a.img.gaussian_beam(err_width, shape=(DIM,DIM))
for bl in field:
    field[bl] *= window
    field[bl] = n.abs(a.img.convolve2d(field[bl], err_kernel))
    field[bl] /= field[bl].max()
field = sum(field.values())
field = field.transpose()
p.imshow(n.log10(field.clip(1e-10)), vmin=-2, vmax=0, extent=(-SIZE/2,SIZE/2,-SIZE/2,SIZE/2), origin='lower')
for bl in data[t]['bls']:
    i,j = a.miriad.bl2ij(bl)
    ax,ay,az = aa.get_baseline(0,i)
    bx,by,bz = aa.get_baseline(i,j)
    print j, ax+bx, ay+by
#    p.plot([ax+bx], [ay+by], '^')
p.show()
