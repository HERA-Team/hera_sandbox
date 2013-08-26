#! /usr/bin/env python

import aipy as a
import numpy as n
import optparse,sys,ephem

o = optparse.OptionParser()
o.set_usage('lst [options] jd1 jd2 ...')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, .1, .1, 1)

print "Filename\t file LST\t aa.sidereal_time(t)"
for F in args:
    uv = a.miriad.UV(F)
    (uvw,t,(i,j)),d = uv.read()
    lst = uv['lst']
    aa.set_jultime(t)
    print F,"\t",ephem.hours(lst),'\t',aa.sidereal_time()


