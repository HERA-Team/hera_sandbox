#!/usr/global/paper/bin/python
"""
input a filename, get out a year-month
"""

import aipy as a, sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('month [options] file1 file2 ...')
o.set_description(__doc__)
opts, args = o.parse_args(sys.argv[1:])

for f in args:
    D=float('.'.join(f.split('.')[1:3]))
    o=a.ephem.Observer()
    o.date = a.phs.juldate2ephem(D)
    print "%02d-%02d"%(o.date.datetime().year-2000,o.date.datetime().month)



