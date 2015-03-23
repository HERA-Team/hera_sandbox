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

npz = n.load('tsys.npz')
tsys = npz['tsys']


sun = a.cal.get_catalog(opts.cal, ['Sun'])['Sun']
lsts = n.arange(0,2*n.pi,LSTRES)
#tsys = n.polyval(tsys_poly, lsts)
#p.plot(lsts,tsys); p.show()
jd_range = n.arange(jd1,jd2,JDSTEP)
print jd1, jd2, jd_range.size

jd_day = 0
lstbin = None

lst_wfall = []
tsys_wfall = []
streak = 0
do_add = False
months = []
for jd in jd_range:
    if int(jd) != jd_day:
        lst_wfall.append(n.zeros_like(lsts))
        jd_day = int(jd)
        #print jd_day
        aa.set_jultime(jd_day)
        tsys_wfall.append(tsys)
        date = str(aa.date).split()[0]
        yr,mo,dy = date.split('/')
        if dy == '1': months.append((mo,yr,jd))
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
        #print jd, streak, float(sun.alt), do_add
    else:
        streak -= 1
        lstbin = (lstbin + 1) % lsts.size
    if do_add: lst_wfall[-1][lstbin] += JDSTEP / a.ephem.second
lst_wfall = n.array(lst_wfall)
tsys_wfall = n.array(tsys_wfall)
print tsys_wfall.max(), tsys_wfall.min()

use = lst_wfall/tsys_wfall**2
use = use.clip(.5 * use.max(),use.max())
#print (n.sqrt(lst_wfall)/tsys_wfall).max()
#print (n.sqrt(lst_wfall)/tsys_wfall).min()
#p.subplot(121)
p.imshow(use, origin='lower', interpolation='nearest', aspect='auto', extent=(0,24,jd1,jd2),cmap='gist_yarg')
#p.imshow(tsys_wfall, origin='lower', interpolation='nearest', aspect='auto', extent=(0,24,jd1,jd2),cmap='hot', alpha=.5)
#p.plot([ 2, 2],[jd1,jd2], 'm--', label=None)
#p.plot([13,13],[jd1,jd2], 'm--', label=None)
for mo,yr,jd in months:
    p.plot([0,24], [jd,jd], label='%s/%s'%(yr,mo))
p.xlim(0,24)
p.xlabel('LST')
p.ylabel('JD')
p.ylim(jd1,jd2)
p.legend(loc='lower center')

#p.subplot(122)
#p.plot(n.linspace(0,24,lsts.size), lst_wfall.sum(axis=0))

p.show()
