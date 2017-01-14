#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys

lstres = 2 * n.pi * 43. / a.const.sidereal_day

aa = a.cal('psa6240_v003', n.array([.150]))
times,dat,flg = {}, {}, {}
lsts = {}
for filename in sys.argv[1:]:
    #times[filename],dat[filename],flg[filename] = C.arp.get_dict_of_uv_data([filename], 'cross', 'xx,yy')
    times[filename],dat[filename],flg[filename] = C.arp.get_dict_of_uv_data([filename], '1_4', 'xx,yy')
    for i,t in enumerate(times[filename]):
        aa.set_jultime(t)
        lst = n.around(aa.sidereal_time() / lstres) * lstres
        lsts[lst] = lsts.get(lst,[]) + (filename,i)
    
lstkeys = lsts.keys()
lstkeys.sort()
for lst in lstkeys:
    plt_x
    print lst, lsts[lst]

