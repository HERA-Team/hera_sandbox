#! /usr/bin/env python
import aipy as a, numpy as n, sys, os, ephem, optparse, glob
from time import strftime, gmtime, mktime, strptime
import pylab as p

def tstr2jd(tstr, ifmt='%m/%d/%y %H:%M:%S', tz=''):
    try:
        #tstr = strftime('%Y/%m/%d %H:%M:%S',
        #    gmtime(mktime(strptime(tstr+' '+tz, ifmt+' %Z'))))
        #tstr = strftime('%Y/%m/%d %H:%M:%S', strptime(tstr, ifmt))
        tstr = strftime('%Y/%m/%d %H:%M:%S', strptime(tstr, ifmt))
        return a.phs.ephem2juldate(ephem.date(tstr))
#    except(ValueError): return []
    except(KeyboardInterrupt): return []

def parse_gom_line(line, filename):
    fields = line.split()
    if len(fields) < 8 and len(fields) > 0:
        ifmt='%m%d%y %H:%M:%S'
        date = tstr2jd(fields[0] + ' ' + fields[1],ifmt = '%Y/%m/%d %H:%M:%S')
        #date = tstr2jd(' '.join([filename[-10:-4], fields[0]]), ifmt=ifmt)
        temps = map(float, fields[2:])
    else:
        date = tstr2jd(' '.join(fields[:2]), ifmt = '%Y/%m/%d %H:%M:%S')
        temps = map(float, fields[2:])
    return [date] + temps

def grid_jd(jds, temps, binsize=10):
    jdbin = binsize * ephem.second
    nbins = int((jds[-1] - jds[0]) / jdbin)
    wgts,bins = n.histogram(jds, bins=nbins)
    dats,bins = n.histogram(jds, weights=temps, bins=nbins)
    return [dats / wgts, bins]



o = optparse.OptionParser()
o.add_option('-t', '--tempdir', dest='tempdir',
    help='Directory containing temperature data from the gainometer.')
opts,args = o.parse_args(sys.argv[1:])

#print opts.tempdir

if opts.tempdir != None:
    data = []
    T1 = []
    T2 = []
    T3 = []
    T4 = []
    dat = []
    files = glob.glob(opts.tempdir + '2011*.txt')
    files.sort() # sorts files in numberical order, i.e. time order.
    for f in files:
        print 'Reading ', f
        lines = [parse_gom_line(L,f) for L in open(f).readlines()] 
        dat.append(L for L in lines if len(L) == 6)
    dat=n.array(dat)
    for d in dat:
        for t in d:
            data.append(t)
    data = n.array(data)
    print data.shape
    tp = data.transpose()
    print tp.shape
    print data.transpose()
    
#    for p in data:
#         T1.append
#        
#            T1.append(grid_jd(t[0],t[1]))
#            T2.append(grid_jd(t[0],t[1]))
#            T3.append(grid_jd(t[0],t[1]))
#            T4.append(grid_jd(t[0],t[1]))
            
#
p.plot(tp[0], tp[1],'k',label='Inside')
p.plot(tp[0], tp[3],'b',label='Outside')
p.title('Temperature vs. Julian Date')
p.xlabel('JD')
p.ylabel('Kelvin')
p.legend()
p.show()



