#! /usr/bin/env python
import numpy as n, pylab as p, aipy as a
import capo as C
LST_COARSEBIN = 1e-2
import sys

g = {}
lsts = []
for f in sys.argv[1:]:
    print f
    _g = n.load(f)
    ants = _g['ants']
    lsts.append(_g['lst'])
    for cnt, i in enumerate(ants):
        for j in ants[cnt+1:]:
            d = []
            for lst in lsts[-1]:
                bij = C.pspec.uv2bin(10*i,10*j,lst, lst_res=LST_COARSEBIN)
                d.append(_g[str(bij)])
            bl = a.miriad.ij2bl(i,j)
            g[bl] = g.get(bl, []) + d

lsts = n.concatenate(lsts)
lsts = n.where(lsts > 3, lsts - 2*n.pi, lsts)
order = n.argsort(lsts)
lsts = lsts[order]
avg = 0
for bl in g:
    g[bl] = n.array(g[bl])[order]
    avg += g[bl]
avg /= len(g)
med = n.array(g.values())
med = n.median(med, axis=0)

p.plot(lsts, avg)
p.plot(lsts, med)
for bl in g:
    norm = g[bl]/avg - 1
    if n.abs(n.average(norm)) > .2:
        print a.miriad.bl2ij(bl), n.average(norm), norm.max(), norm.min()
    p.plot(lsts, g[bl]/avg - 1, '-', label=i, alpha=.2)

p.show()

