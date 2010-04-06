#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import pfits, sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)

times = []
data = {'sum':[]}
for filename in args:
    print filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        d = n.abs(d.take(chans))
        if len(times) == 0 or times[-1] != t:
            times.append(t)
            data['sum'].append(d)
        else: data['sum'][-1] += d
        bl = a.miriad.ij2bl(i,j)
        if not data.has_key(bl): data[bl] = []
        data[bl].append(d)

times = n.array(times)
d = n.array(data['sum'])
for ch in range(d.shape[-1]):
    p.semilogy(times, d[:,ch], label='%d' % (chans[ch]))
del(data['sum'])


#for bl in data:
#    i,j = a.miriad.bl2ij(bl)
#    data[bl] = n.array(data[bl])
#    for ch in range(data[bl].shape[-1]):
#        p.semilogy(data[bl][:,ch], label='%d-%d/%d' % (i,j,chans[ch]))

p.legend()
p.show()
