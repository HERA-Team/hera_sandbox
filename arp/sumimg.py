#! /usr/bin/env python
import aipy as a, sys

mode = 'add'
addargs = []
subargs = []
for arg in sys.argv[1:]:
    if arg == '-a': mode = 'add'
    elif arg == '-s': mode = 'sub'
    elif arg == '-h':
        print 'sumimg.py [-a/-s] file1 ... [-a/-s] more files'
        sys.exit(0)
    else:
        if mode == 'add': addargs.append(arg)
        else: subargs.append(arg)
        
m = None
for img in addargs:
    print 'Adding', img
    d,kwds = a.img.from_fits(img)
    if m is None: m = d
    else: m += d
for img in subargs:
    print 'Subtracting', map
    d,kwds = a.img.from_fits(img)
    if m is None: m = -d
    else: m -= d
a.img.to_fits('sumimg.fits', m, clobber=True, **kwds)
