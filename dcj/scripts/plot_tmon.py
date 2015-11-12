import aipy as a, numpy as n, sys, os, ephem, optparse, glob
import datetime
import time
import matplotlib
from matplotlib import dates
import pylab as pl
from astropy.time import Time

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

times = []
jtimes = []
Tempdata = []
files = args
files.sort() # sorts files in numberical order, i.e. time order.
print "loading files"
for file in files:
    dat = n.loadtxt(file)
    Tempdata.append(dat)
Tempdata = n.vstack(Tempdata)


times = Time(Tempdata[:,0],scale='utc',format='jd')
fig = pl.figure()
ax = pl.subplot(111)
pl.plot_date(times.plot_date,Tempdata[:,1]-273,'b-',label='Container')
pl.plot_date(times.plot_date,Tempdata[:,2]-273,'g-',label='East Balun')
#pl.plot_date(times.plot_date,Tempdata[:,3]-273,'c-',label='East Cable')

pl.xticks(rotation=45)
pl.subplots_adjust(bottom=0.2)
pl.xlabel('UTC')
pl.ylabel('Temperature (C)')
pl.draw()
pl.legend()

# do a nightly average temperature
#hours = times.datetime[3]+2 #add two hours to get local time in SA
nights = list(set(n.round(times.jd)))
avg_nights,balun_lows,balun_means = [],[],[]
t_min = []
hours = n.array([t.hour for t in times.datetime])+2 #local hours
for i,night in enumerate(nights):
    night_time = n.argwhere(n.logical_and(n.round(times.jd)==night,# #select my night
                        n.logical_or(hours>6,
                        hours<18))).squeeze()

    if night_time.size==0:continue
    t_min.append(times[night_time][Tempdata[night_time,2].argmin()].jd)
    avg_nights.append(night)
    balun_lows.append(n.min(Tempdata[night_time,2]))
    balun_means.append(n.mean(Tempdata[night_time,3]))
    print t_min[-1], balun_lows[-1],balun_means[-1]
avg_nights = Time(n.array(avg_nights)+0.5+1/12.,scale='utc',format='jd')
min_times = Time(n.array(t_min),scale='utc',format='jd')
balun_lows = n.array(balun_lows)-273.
balun_means = n.array(balun_means)-273.
night_times = Time(nights,scale='utc',format='jd')
pl.plot_date(min_times.plot_date,balun_lows,label='Balun Low')
#pl.plot_date(avg_nights.plot_date,balun_means,label='Balun Mean')
pl.legend()
pl.show()

