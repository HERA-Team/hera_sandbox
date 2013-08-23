#! /usr/bin/env python

import sys
import waterfall
import datetime
import optparse

o = optparse.OptionParser()
o.add_option('-1','--time1',default=datetime.datetime(2012,11,1,0,0,0),help='start time:  mon/day/year.hr:min')
o.add_option('-2','--time2',default=datetime.datetime(2012,12,1,0,0,0),help='stop time:  mon/day/year.hr:min')
o.add_option('-r','--reftime',default='auto',help="reference time:   mon/day/year.hr:min or 'auto'")
o.add_option('-b','--band',default=(100.0,200.0),help='b1,b2')
o.add_option('-n','--nchan',default=1024,help='n')
o.add_option('-p','--path',default='/Users/daviddeboer1/Documents/Projects/PAPER/data/rfi/RFI_FILES/2012-2013')

opts,args = o.parse_args(sys.argv[1:])

def mktime(ti):
    dt = ti.split('.')
    if len(dt)==1:
        time=('0','0')
    else:
        time = dt[1].split(':')
    date = dt[0].split('/')
    return datetime.datetime(int(date[2]),int(date[0]),int(date[1]),int(time[0]),int(time[1]))


if type(opts.time1) == str:
    opts.time1 = mktime(opts.time1)
if type(opts.time2) == str:
    opts.time2 = mktime(opts.time2)
if type(opts.reftime) == str and opts.reftime!='auto':
    opts.reftime = mktime(opts.reftime)


if type(opts.band) == str:
    ba = opts.band.split(',')
    opts.band = (float(ba[0]),float(ba[1]))


w = waterfall.wf(opts.time1,opts.time2,opts.reftime,opts.nchan,opts.band,opts.path)
w.get_wf()
w.timePlot()
w.scanPlot()
