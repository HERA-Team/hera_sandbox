#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys

aa = a.cal.get_aa(sys.argv[-1], n.array([.150]))
img = a.img.Img(size=300., res=.5)
x,y,z = img.get_top(center=(300,300))
SH = x.shape
x,y,z = x.flatten(), y.flatten(), z.flatten()
bm = aa[0].bm_response((x,y,z), pol='x')**2
bm.shape = SH

p.imshow(bm, interpolation='nearest')
p.show()

