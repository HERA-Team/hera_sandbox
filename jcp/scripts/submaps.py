#!/usr/global/paper/bin/python

"""
Subtracts two Healpix maps from each other.  The input format is fits files.
"""

import aipy as a, ephem as e, numpy as n, sys, optparse, os

o = optparse.OptionParser()
o.set_usage('submap.py [options]')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True)

o.add_option('-m', '--map', dest='map', default='out.fits',
    help='The location to save the output map.')
o.add_option('-i', '--in_map_i', dest='in_map1', 
    help='First map to read in.')
o.add_option('-j', '--in_map_j', dest='in_map2', 
    help='Second map to read in.')
opts,args = o.parse_args(sys.argv[1:])

imap1 = a.map.Map(fromfits=opts.in_map1)
imap2 = a.map.Map(fromfits=opts.in_map2)

h=imap1

h.map.map = (imap1.map.map/imap1.wgt.map)/(imap2.map.map/imap2.wgt.map)
h.wgt.map = n.ones_like(imap1.wgt.map)


print 'Saving to', opts.map
h.to_fits(opts.map)
