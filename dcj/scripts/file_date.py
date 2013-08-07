#!/usr/bin/env python
"""
Read in a list of files and output the date of the file.
"""

import aipy as a,optparse,sys,os,numpy as n
o = optparse.OptionParser()
o.set_usage('file_date.py [options] <files>')
o.set_description(__doc__)
o.add_option('-d',action='store_true',
    help='truncate date to day and print number of files for that day')
o.add_option('-u',action='store_true',
    help="print a list of julian days for all input files")
opts, args = o.parse_args(sys.argv[1:])


O = a.ephem.Observer()
counts = {}
for file in args:
    jd = float('.'.join(os.path.basename(file).split('.')[1:3]))
    O.date = a.phs.juldate2ephem(jd)
    if not (opts.d or opts.u):
        print file, O.date
        continue
    else:
        day = int(jd)
        try:
            counts[day] += 1
        except(KeyError):
            counts[day] = 1
if opts.d:
    for day in counts:
        O.date = a.phs.juldate2ephem(float(day))
        print str(O.date).split()[0],counts[day] #print the date (no time) and file count
if opts.u:
    for day in n.sort(counts.keys()):
        print day
