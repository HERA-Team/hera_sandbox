#! /usr/bin/env python

import sys
import waterfall
import datetime
import argparse

o = argparse.ArgumentParser()
o.add_argument('-1','--time1',default=datetime.datetime(2012,11,1,0,0,0),help='start time:  mon/day/year.hr:min')
o.add_argument('-2','--time2',default=datetime.datetime(2012,11,10,0,0,0),help='stop time:  mon/day/year.hr:min')
o.add_argument('-r','--reftime',default='auto',help="reference time:   mon/day/year.hr:min or 'auto'")
o.add_argument('-b','--band',default=(100.0,200.0),help='b1,b2')
o.add_argument('-n','--nchan',default=1024,help='n')
o.add_argument('-p','--path',default='/Users/daviddeboer1/Documents/Projects/PAPER/data/rfi/RFI_FILES/2012-2013')

args = o.parse_args(sys.argv[1:])

def mktime(ti):
    dt = ti.split('.')
    if len(dt)==1:
        time=('0','0')
    else:
        time = dt[1].split(':')
    date = dt[0].split('/')
    if len(data[2])==2:
        data[2] = '20'+date[2]
    return datetime.datetime(int(date[2]),int(date[0]),int(date[1]),int(time[0]),int(time[1]))


if type(args.time1) == str:
    args.time1 = mktime(args.time1)
if type(args.time2) == str:
    args.time2 = mktime(args.time2)
if type(args.reftime) == str and args.reftime!='auto':
    args.reftime = mktime(args.reftime)


if type(args.band) == str:
    ba = args.band.split(',')
    args.band = (float(ba[0]),float(ba[1]))


w = waterfall.wf(args.time1,args.time2,args.reftime,args.nchan,args.band,args.path)
w.get_wf()
w.timePlot()
w.scanPlot()
