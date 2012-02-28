#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys

m = a.map.Map(fromfits=sys.argv[-1])
i = n.arange(m.wgt.map.size)

wgts = n.where(m.wgt.map == 0, 1, m.wgt.map)
flx = m.map.map / wgts
valid = n.where(m.wgt.map > 1e-2, 1, 0)
flx = flx.compress(valid)
i = i.compress(valid)
strong = n.where(flx > 1, 1, 0)
flx = flx.compress(strong)
i = i.compress(strong)

x,y,z = m.px2crd(i, ncrd=3)
ra, dec = a.coord.eq2radec((x,y,z))

p.plot(ra, dec, '.')
p.xlim(0, 2*n.pi)
p.ylim(-n.pi/2, n.pi/2)
p.xticks(n.arange(0, 2*n.pi, n.pi/12))
p.yticks(n.arange(-n.pi/2, n.pi/2, n.pi/6))
p.grid()
p.show()

