#! /usr/bin/env python
"""
Flag antennae with 2sigma (low) deviantion from mean autocorrelation over time and channel
and
Flag antennae (i) with 1sigma (low) deviation from mean(|V_{ij}|) over time, channel and j
"""
import sys, numpy as np, aipy, optparse, capo
o = optparse.OptionParser()
o.set_usage('meanVij.py [options] *.uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o, cal=False, pol=True)
o.add_option('-c', '--chan', dest='chan', default=None,\
                     help='Channel range in form "lo_hi" to average over.')
o.add_option('-t', '--time', dest='time', default=None,\
                     help='Time range in form "lo_hi" (index) to average over.')
opts,args = o.parse_args(sys.argv[1:])

#get initial info from first file
uv = aipy.miriad.UV(args[0])
nants,nchan = uv['nants'],uv['nchan']
uv.select('antennae',0,1,include=True) #XXX assumes ants 0 and 1 are in uv file
T=[]
for p,d in uv.all(): 
    _,t,_ = p
    T.append(t)
ntimes = len(T)
del(uv);del(T) #scrap potential confusing variables

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
    for i in range(nants):
        for j in range(nants):
            if i==j:
                #get abs(autos)
                try:
                    auto_stor[i,:,:] += np.absolute(data[(i,j)][opts.pol].T)
                    wgt_auto[i,:,:] += np.ones((nchan,ntimes),dtype='complex128')
                except KeyError:
                    print 'KeyError for auto-%i'%i
            else:
                #get abs(crosses)
                try:
                    vis_stor[i,:,:] += np.absolute(data[(i,j)][opts.pol].T)
                    wgt_vis[i,:,:] += np.ones((nchan,ntimes),dtype='complex128')
                except KeyError:
                    if i<j: print 'KeyError on (%i,%i)'%(i,j) #this should not happen
                    continue

#Data analysis
mean_auto = auto_stor/wgt_auto
mean_vis = vis_stor/wgt_vis

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

#array indicies == antenna numbers
print np.where(avg_auto < (mean_avg_auto - 2*std_avg_auto))
print np.where(avg_vis < (mean_avg_vis - std_avg_vis))

