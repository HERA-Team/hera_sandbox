#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p

if True:
    aa = a.cal.get_aa('p32_maxred', n.array([.130,.180]))
    nants = 27
    skip_ants = []
elif True:
    aa = a.cal.get_aa('pgb015_v006', n.array([.130,.180]))
    nants = 16
    skip_ants = [1,7,13,15]
elif True:
    aa = a.cal.get_aa('pgb322_v002_gc', n.array([.130,.180]))
    nants = 32
    skip_ants = []

for i in range(nants):
  for j in range(i+1,nants):
    if i in skip_ants or j in skip_ants: continue
    u,v,w = aa.gen_uvw(i,j,'z')
    u,v = u.flatten(), v.flatten()
    p.plot(u,v,'k-') 
    p.plot(-u,-v,'k-')

rstep = .25
#rs = 10**n.arange(1,2.5,rstep)
rs = 2**(n.arange(3,8,1) +.5)
for r in rs:
    th = n.arange(0, 2*n.pi+.02, .01)
    x,y = r*n.cos(th), r*n.sin(th)
    p.plot(x,y,'r-')

p.grid()
p.xlim(-200,200)
p.ylim(-200,200)
p.xlabel('u',size=14)
p.ylabel('v',size=14)
p.show()
