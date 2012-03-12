#!/usr/bin/env python
import aipy as a, numpy as n, sys
from pylab import *
import optparse


o = optparse.OptionParser()
opts, args = o.parse_args(sys.argv[1:])


k = 160
gather = {}



for filename in args:
    print filename
    uv = a.miriad.UV(filename)
    uv.select('auto',0,0)
    freqs = n.arange(uv['nchan']) * 100.0/256.0 + 100.0
    nants = n.arange(uv['nants']) 
    print ' STUFF:',nants,len(nants),len(freqs)
    cnt = n.zeros( (len(nants),len(freqs)), n.int)
    avg = n.zeros( (len(nants),len(freqs)), n.float)
    max = n.zeros( (len(nants),len(freqs)), n.float)
    min = 100.0 * n.ones( (len(nants),len(freqs)), n.float)
#   r = 0
#   for row in cnt:
#       r += 1
#       print r, row
#   print cnt[0][k]
    tmin = 0
    for p,d,f in uv.all(raw=True):
        crd,t,(i,j) = p
        if not tmin is None: tmin = t
        d = n.abs(d)
#       print k, d[k]
        if i == j and d[k] > 0:
            cnt[i] += 1
            avg[i] += d
            max[i] = n.where(d > max[i], d, max[i])
            min[i] = n.where(d < min[i], d, min[i])
#           print t,max[i][k], min[i][k]
    del(uv)
    gather[t] = {'max':max,'min':min,
        'avg':n.array([avg[i]/cnt[i] for i in range(len(nants))]),
        'cnt':cnt, 'freqs':freqs}
        
ts = n.sort(gather.keys())
#Generate a stackplot of all input files
tmin = ts[0]
tmax = ts[-1]

def jd2hrs(t): #UT hours
    return ((t % int(t))*24+12) % 24
hrs = n.array([jd2hrs(t)+2 for t in ts]) #add 2 tz
dt = n.diff(hrs)[0]
waterfall = []
wf_scale = 10

for t in ts:
    yshift = (1-10**7*(t-tmin)/tmin)
    print t,yshift
    max,min,avg = gather[t]['max'],gather[t]['min'],gather[t]['avg']
    cnt,freqs = gather[t]['cnt'],gather[t]['freqs']
    for i in range(wf_scale): waterfall.append(max[0])
    if jd2hrs(t) % 12 < 1: 
        semilogy(freqs,max[0]*yshift,'k')
        annotate('noon ',(freqs[len(freqs)*2/3],max[0][len(freqs)*2/3]))
    elif jd2hrs(t) % 0 < 1: 
        semilogy(freqs,max[0]*yshift,'k',linewidth=2)
        annotate('midnight '+str(t-2455000),(freqs[len(freqs)*2/3],max[0][len(freqs)*2/3]))
    else: semilogy(freqs,max[0]*yshift,'b',alpha=0.5)
    xlabel('freq MHz')
    ylabel('Power')

waterfall = n.array(waterfall)
figure()
imshow(n.log10(waterfall))
colorbar()
show()
#    for i in range(len(nants)):
#        p.subplot(2, 4, i+1)
#        p.title('ACF %d' %i)
#        p.xlabel('FREQ (MHz)')
#        p.ylabel('POWER/RATIO')
#        p.semilogy(freqs, max[i], 'r-', label='max' )
#        pavg = avg[i] / cnt[i]
#        p.semilogy(freqs, pavg, 'b-', label='avg' )
#        p.semilogy(freqs, min[i], 'b-', label='min' )
#        ratio = (max[i]-pavg) / (pavg-min[i])
#        p.semilogy(freqs, ratio, 'k-', label='ratio')
##       p.legend()
##     print '# type freq average, max, min'
#        p.ylim(0.5,100)
##       for k in range(len(freqs)/10):
##           print i,k, pavg[k], max[i][k], min[i][k], ratio[k]
