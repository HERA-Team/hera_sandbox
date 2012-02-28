#! /usr/bin/env python
import aipy as a, numpy as n, ephem, sys

for f in sys.argv[1:]:
    data, kwds = a.img.from_fits(f)
    ra = ephem.hours(kwds['ra'] / 180. * n.pi)
    dec = ephem.degrees(kwds['dec'] / 180. * n.pi)
    if dec > ephem.degrees('18') \
            and dec < ephem.degrees('58') \
            and ra > ephem.hours('8') \
            and ra < ephem.hours('16'):
        print f, 'RA=%s, DEC=%s' % (ra, dec)
