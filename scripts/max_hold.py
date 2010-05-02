#! /usr/bin/env python
import aipy as a, numpy as n, sys, pylab

data = {}

for filename in sys.argv[1:]:
    print filename
    uv = a.miriad.UV(filename)
    freqs = n.arange(uv['nchan']) * uv['sdf'] + uv['sfreq']
    for p,d,f in uv.all(raw=True):
        crd,t,(i,j) = p
        if i == j: bl = 'auto'
        if i != j: bl = 'cross'
        if not data.has_key(bl):
            data[bl] = {'avg':0, 'max':0, 'clean':0, 'cnt':0, 'clean_cnt':0}
        data[bl]['cnt'] += 1
        data[bl]['clean_cnt'] += n.logical_not(f)
        data[bl]['clean'] += n.abs(n.where(f, 0, d))
        d = n.abs(d)
        data[bl]['avg'] += d
        data[bl]['max'] = n.where(d > data[bl]['max'], d, data[bl]['max'])
    del(uv)

cnt = 0
for bl in data:
    cnt += 1
    pylab.subplot(2, 1, cnt)
    clean = data[bl]['clean'] / data[bl]['clean_cnt']
    avg = data[bl]['avg'] / data[bl]['cnt']
    max = data[bl]['max']
    pylab.semilogy(freqs, max, 'r-', label='%s-max' % bl)
    pylab.semilogy(freqs, avg, 'b-', label='%s-avg' % bl)
    pylab.semilogy(freqs, clean, 'k-', label='%s-cln' % bl)
    pylab.legend()
    print '# type freq clean, average, max'
    for i in range(max.shape[0]):
        print bl, i, clean[i].real, avg[i].real, max[i].real
pylab.show()
