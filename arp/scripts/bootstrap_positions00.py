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
print delays[:2]

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)

peak_delays = {}
srcvec = {}
for src in cat: srcvec[src] = {}

def find_peaks(dly, npeaks=100):
    peaks = {}
    adly = n.abs(dly)
    order = list(n.argsort(adly)[-npeaks:])
    order.reverse()
    for px in order:
        #if peaks.has_key(px-1):
        #    _sum,_wgt = peaks[px-1]
        #    peaks[px-1] = (_sum + delays[px]*adly[px], _wgt + adly[px])
        #elif peaks.has_key(px+1):
        #    _sum,_wgt = peaks[px+1]
        #    peaks[px+1] = (_sum + delays[px]*adly[px], _wgt + adly[px])
        if n.max(adly[max(px-2,0):px+2]) > adly[px]: continue
        else: peaks[px] = (delays[px]*adly[px], adly[px])
    rv = []
    for px in order:
        try:
            _sum,_wgt = peaks[px]
            rv.append(_sum/_wgt)
            #rv.append(_sum/_wgt)
        except(KeyError): pass
    return rv



for filename in args:
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if n.all(f): continue
        if not peak_delays.has_key(t):
            aa.set_jultime(t)
            cat.compute(aa)
            for src in cat: srcvec[src][t] = cat[src].get_crds('top')
        bl = a.miriad.ij2bl(i,j)
        valid = n.logical_not(f).astype(n.float)
        gain = n.sqrt(n.average(valid**2))
        d *= valid
        ker = n.fft.ifft(valid)
        _d = n.fft.ifft(d)
        _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
        _d += info['res'] / gain
        _d = n.abs(_d)
        if not peak_delays.has_key(bl): peak_delays[bl] = {}
        peak_delays[bl][t] = find_peaks(_d)[:len(cat)]
        print t, (i,j), peak_delays[bl][t], [n.dot(srcvec[src][t], aa.get_baseline(i,j)) for src in srclist]

def gen_line(slope, offset):
    dx, dy = offset
    x = n.arange(-2000, 2000, 2)
    y = slope * x
    return x+dx, y+dy

bz_max = 20. # ns
p.plot([0], [0], '.')
for bl in peak_delays:
    i,j = a.miriad.bl2ij(bl)
    ax,ay,az = aa.get_baseline(0,i)
    bx,by,bz = aa.get_baseline(i,j)
    p.plot([ax+bx], [ay+by], '^')
    #print ax+bx,ay+by,(i,j)
    for t in peak_delays[bl]:
        for src,dly in zip(srclist, peak_delays[bl][t]):
            x,y,z = srcvec[src][t]
            dly_err = n.sqrt((delays[1]/2)**2 + (z*bz_max)**2 + 10**2)
            #mag_xy_xyz = n.sqrt(x**2 + y**2) / n.sqrt(x**2 + y**2 + z**2)
            #proj_dly = dly * mag_xy_xyz
            #dx,dy = -x*dly, -y*dly
            #dx,dy = x*dly, y*dly
            #linx,liny = gen_line(y/x, (dx,dy))
            linx,liny = gen_line(-x/y, (0,(dly+dly_err)/y))
            p.plot(linx+ax, liny+ay)
            linx,liny = gen_line(-x/y, (0,(dly-dly_err)/y))
            p.plot(linx+ax, liny+ay)
            #b = dly / (x**2 + y**2)
            #p.plot([-b*x], [-b*y], '.')
            p.plot([1000*x], [1000*y], '.')

p.xlim(-1500,1500)
p.ylim(-1500,1500)
p.show()
