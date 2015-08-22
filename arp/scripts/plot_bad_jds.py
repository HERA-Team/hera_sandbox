#! /usr/bin/env python
import numpy as n, pylab as p
import sys

files = {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    f = n.load(filename)
    for fname,cnt in zip(f['files'],f['cnt']):
        files[fname] = files.get(fname,0) + cnt
    #p.plot(f['jds'], f['bad'], '.')
p.show()

fs = files.keys(); fs.sort()
print fs
#bjd = n.array(jds.keys())
cnt = [files[f] for f in fs]

p.plot(cnt, '.')
p.show()
