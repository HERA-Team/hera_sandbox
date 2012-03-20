import aipy as a, numpy as np, sys, os, ephem, optparse, glob
import datetime
import time
import matplotlib
from matplotlib import dates
matplotlib.use('Agg')
import pylab as pl

arbitrary_offset = 200

def tstr2jd(tstr, ifmt='%Y/%m/%d %H:%M:%S', tz=''):
    tstr = time.strftime('%Y/%m/%d %H:%M:%S', time.strptime(tstr, ifmt))
    return a.phs.ephem2juldate(ephem.date(tstr))

def today():
    cur_time = time.localtime(time.time())
    tstr =  '%04d%02d%02d'%(cur_time.tm_year,cur_time.tm_mon,cur_time.tm_mday)
    return tstr

def format_times(ymd,hms,ifmt ='%Y/%m/%d %H:%M:%S'):
    t = time.strptime(ymd+' '+hms, ifmt)
    return datetime.datetime(t.tm_year,t.tm_mon,t.tm_mday,t.tm_hour,t.tm_min,t.tm_sec)

times = []
jtimes = []
Tin = []
Tout = []

#tempdir = '/home/obs/Temperature_Logs/'
#files = glob.glob(tempdir + today() + '*.txt')
files = sys.argv[1:]
files.sort() # sorts files in numberical order, i.e. time order.
print "loading files"
for file in files:
    dat = np.loadtxt(file,
        dtype={'names':('ymd','hm','t_home','t1','t2'),
        'formats':('S12','S8','float','float','float')})
    for d in dat:
        jtimes.append(tstr2jd(d['ymd']+' '+d['hm']))
        times.append(format_times(d['ymd'],d['hm']))
        Tin.append(d['t_home']-273)
        Tout.append(0.5 * (d['t1']+d['t2'])-273 + arbitrary_offset)

fig = pl.figure()
dstart = '%d/%d/%d'%(times[0].month,times[0].day,times[0].year)
jdstart = np.floor(jtimes[0])
#outpng = '/home/obs/Temperature_Logs/'+str(int(jdstart))+'_daily_temps.png'
outpng = os.path.join(os.path.dirname(sys.argv[-1]),'all_temps.png')
ax = pl.subplot(111)

pl.plot_date(times,Tin,'b-',label='Container')
pl.plot_date(times,Tout,'g-',label='Balun')

Mfmt = dates.DateFormatter('%m - %d')
mfmt = dates.DateFormatter('%H')
ax.xaxis.set_major_locator(dates.DayLocator())
if len(files)<(24*5):
    ax.xaxis.set_minor_locator(dates.HourLocator(interval=int(len(files)*4/(24*3))))
    ax.xaxis.set_minor_formatter(mfmt)

ax.xaxis.set_major_formatter(Mfmt)
#pl.xticks(rotation='vertical')
labels = ax.get_xticklabels() 
for label in labels: 
    label.set_rotation(90) 
pl.title('Temperature Reading since JD %s' % dstart)
pl.xlabel('SAST [Month, Day]')
pl.ylabel('Temperature (C)')
pl.subplots_adjust(bottom=0.2)
pl.draw()
pl.legend(loc='upper left')
fig.savefig(outpng,format='png')
#fig.savefig('/home/obs/Temperature_Logs/today.png')
