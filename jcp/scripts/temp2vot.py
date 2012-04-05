#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys
import atpy, glob, time, datetime, ephem

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

def tstr2jd(tstr, ifmt='%Y/%m/%d %H:%M:%S', tz=''):
    tstr = time.strftime('%Y/%m/%d %H:%M:%S', time.strptime(tstr, ifmt))
    return a.phs.ephem2juldate(ephem.date(tstr))

def format_times(ymd,hms,ifmt ='%Y/%m/%d %H:%M:%S'):
    t = time.strptime(ymd+' '+hms, ifmt)
    return datetime.datetime(t.tm_year,t.tm_mon,t.tm_mday,t.tm_hour,t.tm_min,t.tm_sec)

arbitrary_offset = 200
ofilename = 'xTdata.vot'

t = atpy.Table()

Ttimes = []
jtimes = []
Tin = []
Tout = []
files = args
print 'sorting files'
files.sort() # sorts files in numberical order, i.e. time order.
print 'looping over files'
for filename in files:
    dat = n.loadtxt(filename,
        dtype={'names':('ymd','hm','t_home','t1','t2'),
        'formats':('S12','S8','float','float','float')})
    print 'text loaded'
    for i,d in enumerate(dat):
        if not i % 100:
            jtimes.append(tstr2jd(d['ymd']+' '+d['hm']))
            Ttimes.append(format_times(d['ymd'],d['hm']))
            Tin.append(d['t_home']-273)
            Tout.append(0.5 * (d['t1']+d['t2'])-273 + arbitrary_offset)
t.add_column('Tin', Tin, dtype='<f8',unit='K',description='Inside Labjack')
t.add_column('Tout', Tout, dtype='<f8',unit='K',description='Outside Uncalibrated')
#t.add_column('Ttimes', Ttimes, dtype='<f8',description='Datetime for T')
t.add_column('jtimes', jtimes, dtype='<f8',description='JD for T')

print 'looping over times'
times = n.array(jtimes)
aa = a.cal.get_aa(opts.cal,n.array([.150]))
lst = []
cat = a.src.get_catalog(srcs=['Sun','sgr'])
alt_sun,az_sun = [],[]
alt_sgr,az_sgr = [],[]
for time in times:
    aa.set_jultime(time)
    cat.compute(aa)
    lst.append(aa.sidereal_time() * (12/n.pi))
    alt_sun.append(cat['Sun'].alt * (180/n.pi))
    az_sun.append(cat['Sun'].az * (180/n.pi))
    alt_sgr.append(cat['sgr'].alt * (180/n.pi))
    az_sgr.append(cat['sgr'].az * (180/n.pi))

t.add_column('LST', lst, dtype='<f8', unit='hours')
t.add_column('Alt_Sun', alt_sun, dtype='<f8',unit='deg')
t.add_column('Az_Sun', az_sun, dtype='<f8',unit='deg')
t.add_column('Alt_sgr', alt_sgr, dtype='<f8',unit='deg')
t.add_column('Az_sgr', az_sgr, dtype='<f8',unit='deg')

print 'saving'
t.write(ofilename,overwrite=True)
