#!/bin/sh
''''exec python -u -- "$0" ${1+"$@"} # '''
# vi: syntax=python
import numpy as n, pylab as p
import aipy as a
import sys
import optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])

srclist, cutoff, catalogs = [], None, []
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
aa = a.cal.get_aa(opts.cal, n.array([.15]))


#PLOT_GAIN = False
PLOT_GAIN = 1
colors = 'kkbgcmry'

dly = {}
gain = {}
times = []
lsts = []
for filename in args[1:]:
    print 'Reading', filename
    f = open(filename)
    npz = n.load(f)
    C_phs = npz['C_phs']
    C_amp = 10**npz['C_amp']
    antpos = npz['antpos']
    time = npz['time']
#    print time[0], time[-1]
#    pols = npz['pols']
    pols = ['xx']
    aa.set_jultime(time)
    print aa.sidereal_time()
    lsts.append(aa.sidereal_time()*(12./n.pi))
    times.append(time)
    try:
        print times[-1] - times[-2]
    except: pass
    for pi, pol in enumerate(pols):
        if not dly.has_key(pi): dly[pi],gain[pi] = {},{}
        for i,tau,g in zip(antpos.flatten(), C_phs[pi].flatten(), C_amp[pi].flatten()):
            dly[pi][i] = dly[pi].get(i,[]) + [tau]
            gain[pi][i] = gain[pi].get(i,[]) + [g]
    f.close()

lsts = n.array(lsts)
print lsts
if lsts[0] > lsts[-1]:
    diff = 24 - lsts[0]
    lsts = ((lsts + diff)%24) - diff
print lsts
print times
print antpos
fig = p.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1t= ax1.twiny()
for pi in dly:
#  for i in dly[pi]:
  for ik,i in enumerate(antpos.flatten()):
    c = colors[i/8]
    d = n.array(dly[pi][i])
    g = n.array(gain[pi][i])
    #d_avg = n.average(d)
    d_avg = n.median(d)
    g_avg = n.median(g)
#    print i, 'Dly (avg):', d_avg
    if ik%8==0: print('\n')
#    print '%7.3f'%d_avg,
    print '%7.3f'%g_avg,
    d -= d_avg
    g /= g_avg
#    if len(n.where(g >= 1.4)[0]) > 0: print i,c

    if PLOT_GAIN:
        ax1.plot(times, d, c,label='%d,%s'%(i,pols[pi]))
        ax2.plot(times, g, c, label='%d,%s'%(i,pols[pi]))
        ax1t.plot(lsts, d, c+'.',label='%d,%s'%(i,pols[pi]))
    else:
        p.plot(times, g, c+'.', label='%d,%s'%(i,pols[pi]))

#p.legend()
ax1t.set_xlim(lsts[0],lsts[-1])
#ax1t.xaxis.set_ticks(range(int(lsts[0]),int(lsts[-1])))
ax1.set_xlim(times[0],times[-1])
ax2.set_xlim(times[0],times[-1])
ax1.grid()
ax2.grid()
p.show()
