#! /usr/bin/env python
import aipy as a, sys

mode = 'add'
addargs = []
subargs = []
for arg in sys.argv[1:]:
    if arg == '-a': mode = 'add'
    elif arg == '-s': mode = 'sub'
    elif arg == '-h':
        print 'summap.py [-a/-s] file1 ... [-a/-s] more files'
        sys.exit(0)
    else:
        if mode == 'add': addargs.append(arg)
        else: subargs.append(arg)
        
m = None
for map in addargs:
    print 'Adding', map
    if m is None: m = a.map.Map(fromfits=map)
    else:
        m2 = a.map.Map(fromfits=map)
        m.map.map += m2.map.map
        m.wgt.map += m2.wgt.map
for map in subargs:
    print 'Subtracting', map
    m2 = a.map.Map(fromfits=map)
    m.map.map -= m2.map.map
    m.wgt.map += m2.wgt.map
m.to_fits('summap.fits')
