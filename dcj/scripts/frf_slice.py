#! /usr/bin/env python
import numpy as n, pylab as p,ephem
import sys, aipy as a,re
from pylab import *
import matplotlib.lines as lines
import itertools,optparse
from matplotlib.pyplot import cm
import matplotlib.animation as animation
r2h = 12/n.pi
i,j = 64,49

o = optparse.OptionParser()
o.add_option('--bl',help='Baseline to plot. ex --bl=(1,4)')
o.add_option('--filter_files',type=str,
    help="list of frf filter pkls. 'file1.pkl,file2.pkl'")
opts, args = o.parse_args()
if not opts.filter_files is None:
    filter_files = opts.filter_files.split(',')

i,j = map(int,opts.bl.split(','))
mybl=str(a.miriad.ij2bl(i,j))
print mybl
def file2jd(zenuv):
    return float(re.findall(r'\d+\.\d+', zenuv)[0])
def hms(h):
    return ephem.hours(h)
def dB(x):
    return 10*n.log10(x)
data = {}
time = {} #time is the lst axis... sorry
jd_times = {}
jds = []
for filename in args:
    print filename
    jd = int(n.round(file2jd(filename)))
    t = file2jd(filename)
    try:
        npz = n.load(filename)
    except:
        print "    failed to load"
        continue
    try:
        npz[mybl]

    except(KeyError):
        print "baseline {i}_{j} not found in {f}, skipping".format(i=i,j=j,f=filename)
        continue
    try:
        ntimes = len(npz[mybl])
        print npz['t'+mybl]
        dt = n.diff(npz['t'+mybl])[0]/n.pi #convert radians to fraction of a day
        data[jd].append(npz[mybl])
        time[jd].append(npz['t'+mybl])
        jd_times[jd].append(n.arange(0,ntimes)*dt + file2jd(filename))

    except(KeyError):
	try:
            data[jd] = [npz[mybl]]
            time[jd] = [npz['t'+mybl]]
            jd_times[jd] = [n.arange(0,ntimes)*dt + file2jd(filename)]
        except(KeyError):
            continue
if len(data)==0:
    print "ERROR: no data found"
    sys.exit()
for jd in data:
    data[jd] = n.concatenate(data[jd])#turna list of arrays into a single array
    time[jd] = n.concatenate(time[jd])
    jd_times[jd] = n.concatenate(jd_times[jd])
jd_samples = n.concatenate([jd_times[jd] for jd in jd_times])
for jd in data:
    p.plot(time[jd],n.abs(data[jd]),'k')



p.show()

