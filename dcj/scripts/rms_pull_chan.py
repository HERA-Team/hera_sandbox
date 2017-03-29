import numpy as n
import sys,aipy as a,re
from pylab import *
i,j = 64,49
mybl=str(a.miriad.ij2bl(i,j))
def file2jd(zenuv):
    return float(re.findall(r'\d+\.\d+', zenuv)[0])
times = []
rmss = []
maxes = []
for arg in sys.argv[1:]:
    npz = n.load(arg)
    try:
        rmss.append(n.std(npz[mybl]))
        maxes.append(n.max(npz[mybl]))
    except(KeyError):
        continue
    times.append(file2jd(arg))
    if rmss[-1]>0.03: print arg,rmss[-1]
times = n.array(times)
plot(times-times.min(),rmss)
#plot(times-times.min(),maxes)
text(0.92,-0.07,"+%i"%times.min(),fontsize=10,transform=gca().transAxes)
show()
