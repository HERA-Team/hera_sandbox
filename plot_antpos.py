#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

th = n.arange(0, 2*n.pi, .01)

aa = a.cal.get_aa(opts.cal, .1, .1, 1)
antpos = [aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))]
antpos = n.array(antpos) * a.const.len_ns / 100.
x,y,z = antpos[:,0], antpos[:,1], antpos[:,2]
x -= n.average(x)
y -= n.average(y)
p.plot(x,y, 'k.')
for ant,(xa,ya,za) in enumerate(zip(x,y,z)):
    hx,hy = 10*za*n.cos(th)+xa, 10*za*n.sin(th)+ya
    if za > 0: fmt = '#eeeeee'
    else: fmt = '#a0a0a0'
    p.fill(hx,hy, fmt)
    p.text(xa,ya, str(ant))
p.grid()
#p.xlim(-150,150)
p.xlabel("East-West Antenna Position (m)")
p.ylabel("North-South Antenna Position (m)")
#p.ylim(-150,150)
a = p.gca()
a.set_aspect('equal')
p.show()
