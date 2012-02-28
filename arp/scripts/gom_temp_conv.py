#!/usr/bin/env python
import sys
from time import strftime, gmtime, mktime, strptime

def tstr2utc(tstr, ifmt='%m/%d/%y %H:%M:%S', tz='EDT'):
    return strftime('%m/%d/%y %H:%M:%S',
            gmtime(mktime(strptime(tstr+' '+tz, ifmt+' %Z'))))

dat = {}
for f in sys.argv[1:]:
    print 'Reading', f
    for L in open(f).readlines():
        L = L.split()
        if len(L) != 6: continue
        L[:2] = tstr2utc(' '.join(L[:2])).split()
        if not dat.has_key(L[0]): dat[L[0]] = []
        dat[L[0]].append('\t'.join(L))

for k in dat:
    filename = 'Temp'+k.replace('/','')+'.txt'
    print 'Writing', filename
    f = open(filename, 'w')
    f.write('\n'.join(dat[k]))
    f.close()
