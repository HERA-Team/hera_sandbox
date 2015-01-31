#! /usr/bin/env python
import numpy as n, pylab as p
import sys, aipy as a,re
from pylab import *
r2h = 12/n.pi
i,j = 64,49
mybl=str(a.miriad.ij2bl(i,j))
jdjump = 2456949.
print mybl
def file2jd(zenuv):
    return float(re.findall(r'\d+\.\d+', zenuv)[0])
data = {}
time = {}
jds = []
for filename in sys.argv[1:]:
    jd = int(n.round(file2jd(filename)))
    print 'Reading', filename
    try:
        npz = n.load(filename)
    except:
        print "    failed to load"
    try:
        data[jd].append(npz[mybl])
        time[jd].append(npz['t'+mybl])
    except(KeyError):
	try:
            data[jd] = [npz[mybl]]
            time[jd] = [npz['t'+mybl]]
        except(KeyError):
            continue
for jd in data:
    data[jd] = n.concatenate(data[jd])#turna list of arrays into a single array
    time[jd] = n.concatenate(time[jd])

#make my lst grid    
lsts = n.concatenate([time[jd] for jd in time])
lstmin,lstmax,dlst = lsts.min(),lsts.max(),n.min(n.abs(n.diff(lsts)))
lst_grid = n.arange(lstmin*0.95,lstmax*1.05,dlst*0.8) #note I am tweaking the grid size a little bit after looking at the lst x day plot below

#make a jd grid
jds = n.array(time.keys())
jd_grid = n.arange(jds.min()-1,jds.max()+2,1)
print jd_grid.min(),jds.min(),jds.max(),jd_grid.max()
#grid the data into lst,jd bins
gridded_data = n.zeros((lst_grid.size,jd_grid.size)).astype(n.complex64) #dimensions lstbins x jds
counts = n.zeros_like(gridded_data).astype(int)
jd_inds = n.digitize(jds,jd_grid)
for i,jd in enumerate(jds):
    lst_inds = n.digitize(time[jd],lst_grid)
    gridded_data[lst_inds,jd_inds[i]] += data[jd]
    counts[lst_inds,jd_inds[i]] += n.ones_like(data[jd])
gridded_data[counts>0] /= counts[counts>0]
jds = n.array(time.keys())
imshow(n.abs(gridded_data),aspect='auto',vmax=0.05,interpolation='nearest',extent=(0,jds.max()-jds.min(),lst_grid.max()*r2h,lst_grid.min()*r2h))
colorbar()
text(0.92,-0.07,"+%i"%jds.min(),fontsize=10,transform=gca().transAxes)
#hist(n.abs(gridded_data.ravel()),bins=100)
figure()
imshow(counts,aspect='auto',cmap='Greys',interpolation='nearest',extent=(0,jds.max()-jds.min(),lst_grid.max()*r2h,lst_grid.min()*r2h))
text(0.92,-0.07,"+%i"%jds.min(),fontsize=10,transform=gca().transAxes)

colorbar()
show()
