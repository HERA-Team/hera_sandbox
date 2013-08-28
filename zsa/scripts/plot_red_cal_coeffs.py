#!/bin/sh
''''exec python -u -- "$0" ${1+"$@"} # '''
# vi: syntax=python
import numpy as n, pylab as p
import sys

#PLOT_GAIN = False
PLOT_GAIN = 0
colors = 'kkbgcmry'

dly = {}
gain = {}
times = []
for filename in sys.argv[1:]:
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


print antpos
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
        p.subplot(211); p.plot(times, d, c+',',label='%d,%s'%(i,pols[pi]))
        p.subplot(212); p.plot(times, g, c, label='%d,%s'%(i,pols[pi]))
    else:
        p.plot(times, g, c+'.', label='%d,%s'%(i,pols[pi]))

#p.legend()
p.show()
