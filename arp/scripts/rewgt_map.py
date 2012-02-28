#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

for map in args:
    print 'Reading', map
    m = a.map.Map(fromfits=map)
    wgt_clip = n.where(m.wgt.map == 0, 1, m.wgt.map)
    m.map.map /= wgt_clip
    m.wgt.map /= wgt_clip
m.to_fits('rewgt_'+map, clobber=True)
