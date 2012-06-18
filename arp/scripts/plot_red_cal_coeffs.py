#! /usr/bin/env python
import numpy as n, pylab as p
import sys

PLOT_GAIN = False
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
    pols = npz['pols']
    times.append(time)
    for pi, pol in enumerate(pols):
        if not dly.has_key(pi): dly[pi],gain[pi] = {},{}
        for i,tau,g in zip(antpos.flatten(), C_phs[pi].flatten(), C_amp[pi].flatten()):
            dly[pi][i] = dly[pi].get(i,[]) + [tau]
            gain[pi][i] = gain[pi].get(i,[]) + [g]
    f.close()

for pi in dly:
  for i in dly[pi]:
    c = colors[i/4]
    d = n.array(dly[pi][i])
    #d_avg = n.average(d)
    d_avg = n.median(d)
    print i, 'Dly (avg):', d_avg
    d -= d_avg
    
    if PLOT_GAIN:
        p.subplot(211); p.plot(times, d, c, label='%d,%s'%(i,pols[pi]))
        p.subplot(212); p.plot(times, gain[i], c, label='%d,%s'%(i,pols[pi]))
    else:
        p.plot(times, d, c, label='%d,%s'%(i,pols[pi]))

p.legend()
p.show()
