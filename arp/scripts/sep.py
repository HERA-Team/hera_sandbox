#! /usr/bin/env python

import numpy as n, aipy as a
import sys

aa = a.cal.get_aa(sys.argv[-1], n.array([.150]))
layout = aa.ant_layout

bad = [5, 7, 8, 15, 16, 17, 24, 26, 27, 28, 29, 37, 38, 46, 48, 50, 51, 53, 55, 63, 68, 69, 72, 74, 76, 77, 82, 83, 84, 85, 92, 107, 110]
dr,dc = map(int, sys.argv[-2].split(','))
rv = []
for x in range(layout.shape[0]):
    for y in range(layout.shape[1]):
        i = layout[x,y]
        if x+dr < 0 or y+dc < 0: continue
        try: j = layout[x+dr,y+dc]
        except(IndexError): continue
        if not i in bad and not j in bad: continue
        rv.append('%d_%d' % (i,j))
#print ','.join(rv)
print ' '.join(rv)
