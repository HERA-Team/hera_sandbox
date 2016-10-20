#! /usr/bin/env python
"""
Flag antennae with 2sigma (low) deviantion from mean autocorrelation over time and channel
and
Flag antennae (i) with 1sigma (low) deviation from mean(|V_{ij}|) over time, channel and j
"""
import sys, numpy as np, aipy, optparse, capo, re

o = optparse.OptionParser()
o.set_usage('meanVij.py [options] *.uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o, cal=False, pol=True)
o.add_option('-c', '--chan', dest='chan', default=None,\
                     help='Channel range in form "lo_hi" to average over.')
o.add_option('-t', '--time', dest='time', default=None,\
                     help='Time range in form "lo_hi" (index) to average over.')
o.add_option('--ex_ants',dest='ex_ants',help='Comma-separated list of antennae to exclude.')
o.add_option('-v', '--verbose', dest='verb', action='store_true', help='Toggle verbosity.')
o.add_option('--plot', dest='plot', action='store_true', help='Toggle plotting.')
opts,args = o.parse_args(sys.argv[1:])

def file2jd(zenuv): return re.findall(r'\d+\.\d+', zenuv)[0]

#get initial info from first file
uv = aipy.miriad.UV(args[0])
djd = file2jd(args[0])
JD = int(djd.split('.')[0])
nants,nchan = uv['nants'],uv['nchan']
uv.select('antennae',0,1,include=True) #XXX assumes ants 0 and 1 are in uv file
T=[]
for p,d in uv.all(): 
    _,t,_ = p
    T.append(t)
ntimes = len(T)
del(uv);del(T) #scrap potential confusing variables

ex_ants = map(int,opts.ex_ants.split(','))
valid_ants = [x for x in xrange(nants) if x not in ex_ants]

#parse options
if not opts.time is None: tlo,thi = map(int,opts.time.split('_'))
else: tlo,thi = 0,ntimes-1
if not opts.chan is None: clo,chi = map(int,opts.chan.split('_'))
else: clo,chi = 0,nchan-1

#we stack data in these arrays
vis_stor = np.zeros((nants,nchan,ntimes),dtype='complex128')
auto_stor = np.zeros((nants,nchan,ntimes),dtype='complex128')
wgt_vis,wgt_auto = np.zeros_like(vis_stor),np.zeros_like(auto_stor)

#Data collection
#stack data per uv file to avoid MemoryError(s)
for uv in args:
    print '    Reading %s'%uv
    times,data,flags = capo.arp.get_dict_of_uv_data([uv],'all',opts.pol)
    for i in valid_ants:
        for j in valid_ants:
            if i==j:
                #get abs(autos)
                try:
                    auto_stor[i,:,:] += np.absolute(data[(i,j)][opts.pol].T)
                    wgt_auto[i,:,:] += np.ones((nchan,ntimes),dtype='complex128')
                except KeyError:
                    if opts.verb: print 'KeyError for auto-%i'%i
            else:
                #get abs(crosses)
                try:
                    vis_stor[i,:,:] += np.absolute(data[(i,j)][opts.pol].T)
                    wgt_vis[i,:,:] += np.ones((nchan,ntimes),dtype='complex128')
                except KeyError:
                    if i<j: 
                        if opts.verb:
                            print 'KeyError on (%i,%i)'%(i,j) #this should not happen
                    continue


#Data analysis
mean_auto = auto_stor/wgt_auto
mean_vis = vis_stor/wgt_vis
mean_vis[vis_stor==0.+0.j] = np.nan
mean_auto[vis_stor==0.+0.j] = np.nan

if opts.time is None and opts.chan is None:
    avg_auto = np.nanmean(mean_auto,axis=(1,2))
    avg_vis = np.nanmean(mean_vis,axis=(1,2))
else:
    #there must be a better way to do this
    avg_auto,avg_vis = np.zeros((nants)),np.zeros((nants))
    for i in range(nants):
        avg_auto[i] = np.nanmean(mean_auto[i,clo:chi,tlo:thi])
        avg_vis[i] = np.nanmean(mean_vis[i,clo:chi,tlo:thi])

#terrible variable names, but that's what they are!
mean_avg_auto, mean_avg_vis = np.nanmean(avg_auto),np.nanmean(avg_vis)
std_avg_auto, std_avg_vis = np.nanstd(avg_auto),np.nanstd(avg_vis)
if opts.verb:
    print mean_avg_auto,std_avg_auto
    print mean_avg_vis,std_avg_vis
#array indicies != antenna numbers when ex_ants is specified
ba = np.where(avg_auto < (mean_avg_auto - 2*std_avg_auto))[0]
bv = np.where(avg_vis < (mean_avg_vis - std_avg_vis))[0]
badants = np.concatenate((ba,bv))
badants = set(badants)

print JD, sorted(list(badants))

if not opts.plot: sys.exit(0)

from matplotlib import pyplot as plt
plt.plot(range(nants),avg_vis,'bo')
plt.plot(ex_ants,np.zeros_like(ex_ants),'rs')
plt.axhline(mean_avg_vis)
plt.fill_between(range(nants),mean_avg_vis+std_avg_vis,mean_avg_vis-std_avg_vis,alpha=0.5)
plt.show()

plt.plot(valid_ants,avg_vis[valid_ants]/avg_auto[valid_ants],'bo')
for i in valid_ants:
    if i in badants:
        plt.plot(i,avg_vis[i]/avg_auto[i],'kx',ms=15)
plt.grid()
plt.show()
