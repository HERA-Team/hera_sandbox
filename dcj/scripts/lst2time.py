#! /usr/bin/env python
"""
Recompute the times in a miriad file given the listed sidereal times.
usage lst2tim.py z*uv
"""

import aipy as a
import numpy as n
import optparse,sys,ephem,os

o = optparse.OptionParser()
o.set_usage('lst [options] jd1 jd2 ...')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True)
o.add_option('--lst_res',default=90.,type=n.float,
    help="allowable difference between lst and that computed from listed time [m]")
opts, args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, .1, .1, 1)

def lst2date(lst,date):
    "Output the next julian date corresponding to input lst forward from input date."
    obj = a.fit.RadioFixedBody(lst,0)
    obj.compute(aa)
    aa.set_jultime(n.floor(date)) #start the day at noon UTC
    aa.update()
    aa.date =  aa.next_transit(obj)
    aa.update()
    return aa.get_jultime()
def mfunc(uv,preamble,data,flags):
    uvw,t,(i,j) = preamble
    t = lst2date(uv['lst']-opts.lst_res * 2 * n.pi/a.const.sidereal_day,2456418)
#    t = lst2date(uv['lst'],2456418)
    return (uvw,t,(i,j)),data,flags

for F in args:
    uv = a.miriad.UV(F)
    (uvw,t,(i,j)),d = uv.read()
    lst = uv['lst']
    aa.set_jultime(t)
    print F,
    if n.abs(lst-aa.sidereal_time())>(opts.lst_res*2*n.pi/a.const.sidereal_day):
        print "applying LST correction of:",
        print ephem.hours(lst-aa.sidereal_time()-opts.lst_res * 2 * n.pi/a.const.sidereal_day)
        print ephem.hours(lst-aa.sidereal_time())
        outfile = os.path.basename(F)
        uvi = a.miriad.UV(F)
        print "\t",outfile,
        if os.path.exists(outfile):
            print "exists"
            continue
        else:
            print
        uvo = a.miriad.UV(outfile,status="new")
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi,mfunc=mfunc,raw=True,
            append2hist="lst2time: correction = %s"%(str(ephem.hours(lst- aa.sidereal_time())))
            )

