#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True, max=True, drng=True)
opts, args = o.parse_args(sys.argv[1:])

bm_poly = [ -2.02131219,  11.89480783]
pb_poly = [1.02854332e+09, -9.49707493e+08, 3.64775002e+08,
    -7.46038156e+07, 8.56951433e+06, -5.24246222e+05, 1.33464786e+04]

def ch2freq(chan, sfreq=0.1212890625, sdf=0.00029296875):
    return (sfreq + chan*sdf) * 1e9

def T(Jy, chan=162.5, sfreq=0.1212890625, sdf=0.00029296875):
    freq = ch2freq(chan, sfreq=sfreq, sdf=sdf)
    lam = a.const.c / freq
    #bm = 10**n.polyval(bm_poly, n.log10(freq))
    bm = 10**(n.polyval(pb_poly, n.log10(freq))).clip(-1,0)
    #bm = .4
    print bm
    return Jy * 1e-23 / bm * lam**2 / (2 * a.const.k)

uv = a.miriad.UV(args[0])
sdf, sfreq, nchan, inttime = uv['sdf'], uv['sfreq'], uv['nchan'], uv['inttime']
chans = n.arange(nchan)
if not opts.cal is None:
    aa = a.cal.get_aa(opts.cal, sdf, sfreq, nchan)
else:
    print 'No cal file provided: using data scaling'
    aa = None
del(uv)

data = {}
curtime = None
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)),d in uv.all():
        bl = a.miriad.ij2bl(i,j)
        #if t != curtime:
        #    curtime = t
        #    if not aa is None: aa.set_jultime(t)
        #if not aa is None:
        #    d = T(d / aa.passband(i,j), chans, sfreq, sdf)
        #d = n.abs(d)
        try: data[bl].append(d)
        except(KeyError): data[bl] = [d]

for bl in data: data[bl] = n.ma.array(data[bl])

for cnt,bl in enumerate(data):
    p.subplot(len(data), 1, cnt+1)
    i,j = a.miriad.bl2ij(bl)
    d = data[bl]
    dd = (d[1:-1] - .5 * (d[:-2] + d[2:])) / n.sqrt(3./2)
    dd *= n.sqrt(2 * sdf * 1e9 * inttime)
    if aa is None:
        if i != j:
            bli = a.miriad.ij2bl(i,i)
            blj = a.miriad.ij2bl(j,j)
            d = n.sqrt(data[bli] * data[blj])
        dd /= d[1:-1]
    else:
        dd = T(dd / aa.passband(i,j), chans, sfreq, sdf)
    
    if True:
        dd = n.log10(n.abs(dd))
        if opts.max is None: mx = dd.max()
        else: mx = opts.max
        if opts.drng is None: mn = max(dd.min(), -10)
        else: mn = mx - opts.drng
        p.imshow(dd, vmax=mx, vmin=mn, aspect='auto')
        p.colorbar(shrink=.5)
    elif False:
        h, bins = n.histogram(d[:,600:620], bins=80, range=(-20,20))
        p.plot(.5*(bins[:-1]+bins[1:]), h)
    else:
        p.semilogy(n.ma.std(dd, axis=0))
p.show()
