#! /usr/bin/env python
import numpy as n, pylab as p
import sys

g = {}
lsts = []
for f in sys.argv[1:]:
    print f
    _g = n.load(f)
    for i in _g.files:
        try: g[int(i)] = g.get(int(i),[]) + [_g[i]]
        except(ValueError):
            print f, _g[i]
            lsts.append(_g[i])

lsts = n.concatenate(lsts)
lsts = n.where(lsts > 3, lsts - 2*n.pi, lsts)
order = n.argsort(lsts)
lsts = lsts[order]
avg = 0
for i in g:
    g[i] = n.concatenate(g[i])[order]
    avg += g[i]
avg /= len(g)

p.plot(lsts, avg)
for i in g:
    norm = g[i]/avg - 1
    print i, n.average(norm), norm.max(), norm.min()
    p.plot(lsts, g[i]/avg - 1, label=i, alpha=.5)

p.show()

