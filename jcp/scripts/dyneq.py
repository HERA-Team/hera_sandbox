#! /usr/bin/env python
import numpy as np, aipy as a, os, time, math

aa = a .phs.ArrayLocation(('38:25:59.24',  '-79:51:02.1')) # GB
#aa = a.phs.ArrayLocation(('-30:43:17.4', '21:25:41.9')) # SA

#time_table = np.arange(24)
eq_table = np.arange(24)

syst = time.time()
jd = (syst/86400.) + 2440587.5
aa.set_jultime(jd)
lst = aa.sidereal_time()
hour = 24 * (np.array(lst)/(2 * math.pi))
hourint = int(hour - hour % 1)
if hour % 1 > 0.5: hourint+=1
eq = eq_table[hourint]

os.system("testbed3.py -i %d" % eq)
print "its happened!!"
