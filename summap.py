#! /usr/bin/env python
import aipy as a, sys

m = None
for map in sys.argv[1:]:
    print 'Reading', map
    if m is None: m = a.map.Map(fromfits=map)
    else:
        m2 = a.map.Map(fromfits=map)
        m.map.map += m2.map.map
        m.wgt.map += m2.wgt.map
m.to_fits('summap.fits')
