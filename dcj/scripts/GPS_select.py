#! /usr/bin/env python
import optparse,os,sys
from astropy import time #note: as of 22 July 2013 requires git version of astropy
import numpy as n
o = optparse.OptionParser()
o.add_option('-d',dest='date',type='str',
    help="File selection date [format=iso, eg 2013-01-30, default=today (UTC)].")
opts, args = o.parse_args(sys.argv[1:])
if opts.date is None:
    date = time.Time.now()
else:
    date=time.Time(opts.date,format='iso',scale='utc')
for F in args:
    try:
        t = time.Time(float(F[:-4]),format='gps',scale='utc')
    except(ValueError):continue
    if n.abs((t-date).jd)<0.5: #select files within 12 hours of midnight UTC
        print F

