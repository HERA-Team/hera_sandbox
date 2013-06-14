#!/usr/bin/env python

import urllib as U, numpy as n, datetime as D
from pylab import *
#from collections import Counter
import time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.mlab as mlab
import matplotlib.cbook as cbook

yearsloc    = mdates.YearLocator()   # every year
monthsloc   = mdates.MonthLocator()  # every month
daysloc     = mdates.WeekdayLocator(byweekday=mdates.MONDAY)
yearsFmt = mdates.DateFormatter('%Y')
daysFmt = mdates.DateFormatter('%a %d')
wdloc = mdates.WeekdayLocator(byweekday=mdates.FRIDAY)
wdFmt = mdates.DateFormatter('%a %b %d ')

F = open(sys.argv[-1])
lines = F.readlines()[1:]
lines = [line.split() for line in lines]
dates = {}
for line in lines:
    if line[0].startswith('1969'):continue
    thisday = time.strptime(line[0],"%Y/%m/%d")
    thistime  = datetime.date(thisday.tm_year,thisday.tm_mon,thisday.tm_mday)
    try:dates[thistime] += float(line[1])
    except(KeyError): dates[thistime] = float(line[1])
fig = figure()
ax = fig.add_subplot(111)

days = sort(dates.keys())
files = [dates[day] for day in days]
stem(days,files)
ylabel('files/day')
#print files.shape
fig.autofmt_xdate(bottom=0.2,rotation=45)
show()
sys.exit()


