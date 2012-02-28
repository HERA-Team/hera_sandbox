#! /usr/bin/env python
import aipy as a, numpy as n, sys

for filename in sys.argv[1:]:
    print 'Flattening', filename
    m = a.map.Map(fromfits=filename)
    m.reset_wgt()
    m.map.map = n.abs(m.map.map)
    m.to_fits('flat_' + filename)
