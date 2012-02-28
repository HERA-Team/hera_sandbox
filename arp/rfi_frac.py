#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys

cnt, flag = 0, 0
freqs = None
color = (.8,.8,.8)
for filename in sys.argv[1:]:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    if freqs is None:
        freqs = n.arange(uv['nchan'], dtype=n.float) * uv['sdf'] + uv['sfreq']
    uv.select('antennae', 0, 1)
    for stuff,d,f in uv.all(raw=True):
        cnt += 1
        flag += f

bin_cnt = 1 - flag.astype(n.float)/flag.max()
p.fill_between(1e3*freqs, bin_cnt, edgecolor='black', facecolor=color)
p.xlim(110,180)
p.ylim(.6,1.1)
p.show()
