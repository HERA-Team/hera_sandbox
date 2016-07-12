#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n
import optparse, sys

LSTRES = .01
JDSTEP = LSTRES * a.const.sidereal_day / (2*n.pi) * a.ephem.second

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))
try:
    jd1,jd2 = map(float, args)
except(ValueError):
    jd1,jd2 = map(lambda x: a.phs.ephem2juldate(a.ephem.date(x)), args)


sun = a.cal.get_catalog(opts.cal, ['Sun'])['Sun']
lsts = n.arange(0,2*n.pi,LSTRES)
jd_range = n.arange(jd1,jd2,JDSTEP)
print jd1, jd2, jd_range.size

jd_day = 0
lstbin = None
lst_wfall = []
streak = 0
do_add = False
for jd in jd_range:
    if int(jd) != jd_day:
        lst_wfall.append(n.zeros_like(lsts))
        jd_day = int(jd)
        #print jd_day
    if streak == 0:
        aa.set_jultime(jd)
        sun.compute(aa)
        do_add = (sun.alt < 0)
        new_lstbin = int(n.around(aa.sidereal_time()/LSTRES))
        if not lstbin is None and new_lstbin - lstbin > 1 and do_add:
            lst_wfall[-1][lstbin+1] += JDSTEP / a.ephem.second
        lstbin = new_lstbin    
        if n.abs(sun.alt) > 0.01:
            streak = sun.rise_time - aa.date
            if streak < 0: streak = sun.set_time - aa.date
            if streak < 0: streak = sun.rise_time - aa.date + 1
            streak = int(streak / JDSTEP * .9)
        print jd, streak, float(sun.alt), do_add
    else:
        streak -= 1
        lstbin = (lstbin + 1) % lsts.size
    if do_add: lst_wfall[-1][lstbin] += JDSTEP / a.ephem.second
lst_wfall = n.array(lst_wfall)

p.subplot(121)
p.imshow(lst_wfall, origin='lower', interpolation='nearest', aspect='auto', extent=(0,24,jd1,jd2))

p.subplot(122)
p.plot(n.linspace(0,24,lsts.size), lst_wfall.sum(axis=0))

p.show()
