#! /usr/bin/env python

import numpy as n, aipy as a
import sys

aa = a.cal.get_aa(sys.argv[-1], n.array([.150]))
layout = aa.ant_layout

dr,dc = map(int, sys.argv[-2].split(','))
rv = []
for x in range(layout.shape[0]):
    for y in range(layout.shape[1]):
        i = layout[x,y]
        if x+dr < 0 or y+dc < 0: continue
        try: j = layout[x+dr,y+dc]
        except(IndexError): continue
        rv.append('%d_%d' % (i,j))
print ','.join(rv)
