#! /usr/bin/env python
"""
compute abs(night1/night2) +/- std
averaged over times and channels
"""
import aipy as a, numpy as n, os, sys, glob
import optparse
from pylab import *
o = optparse.OptionParser()
a.scripting.add_standard_options(o,chan=True,cal=True)
o.add_option('-b', '--bin', dest='bin', type='float', default=0.001,
    help='Bin size in LST.  Default 0.001')
o.add_option('-t', '--tempdir', dest='tempdir',
    help='Directory containing temperature data from the (labjack) gainometer.')
o.add_option('--flag_ant',default='',help='Antennae to flag')
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal,n.array([0.15]))


if len(opts.flag_ant):
    flag_ant = ',-'+',-'.join(opts.flag_ant.split(','))
else:
    flag_ant=  ''
def f2jd(f):
    return float('.'.join(f.split('.')[1:3]))
def nighti(f,jds):
    for i,t in enumerate(jds):
        if f.find(str(int(t)))>0:
            return i
    return None
lsts = n.arange(0, 2*n.pi, opts.bin)
print lsts.size
#get the two night jds
nights = list(set([int(f2jd(os.path.basename(f))) for f in args]))
uv = a.miriad.UV(args[0])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
nants = uv['nants']
nchan = len(chans)
del(uv)
data = n.ma.zeros((nants,lsts.size,len(nights),len(chans))).astype(n.complex64)
print data.shape
print "working on ",nights
for filename in args:
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv,'autos'+flag_ant,'xx,yy')
    print filename
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        lbin = int(uv['lst'] / opts.bin)
        lst = lbin * opts.bin
        night = nighti(filename,nights)
        data[i,lbin,night,:] = n.ma.masked_where(f[chans],d[chans])
data = n.ma.masked_where(data==0,data)
R = n.ma.masked_invalid(n.abs(data[:,:,0,:]/data[:,:,1,:]))
GLOBAL_G = n.ma.mean(n.ma.mean(R,axis=0),axis=1)#avg over ants and chans
print GLOBAL_G.shape,lsts.shape

print "amps = {"
for i in range(nants):
    try:
        correction = n.ma.mean(R[i,:,:])
        if n.isnan(1/correction):
            raise(AttributeError)
        newamp = aa[i].amp *n.sqrt(correction)
        print '%d:%7.5f, # (calibrated to %d, a %4.2f%% correction)'%(i,newamp,nights[1],
                100*(newamp - aa[i].amp)/aa[i].amp)
    except(AttributeError):
        print '%d:7.5f, # (not calibrated, listed as flagged by compare_nights.py)'%(i,aa[i].amp)
print "}"
sys.exit()
figure(1)
plot(lsts*12/n.pi,GLOBAL_G)
#GLOBAL_G = n.tile(GLOBAL_G,(1,nchan))#extend back over chans
GLOBAL_G = n.array([GLOBAL_G]*nchan).T
print GLOBAL_G.shape
print R.max(),R.min()
figure()
for i in range(data.shape[0]):
    subplot(8,8,i+1)
    errorbar(chans,n.ma.mean(R[i,:,:]-GLOBAL_G,axis=0),yerr=n.ma.std(R[i,:,:]-GLOBAL_G,axis=0))
figure()
for i in range(data.shape[0]):
    subplot(8,8,i+1)
    errorbar(lsts*12/n.pi,n.ma.mean(R[i,:,:]-GLOBAL_G,axis=1),yerr=n.ma.std(R[i,:,:]-GLOBAL_G,axis=1))
for i in range(data.shape[0]):
    try:
#        subplot(8,8,i+1)
#        plot(lsts*12/n.pi,R[i,:,:],'.')
        print i,n.round(n.ma.mean(R[i,:,:]),2),'+/-',n.round(n.ma.std(R[i,:,:]),4)
    except(AttributeError):
        print '!'
        continue
print "overall"
print n.round(n.ma.mean(R),2),'+/-',n.round(n.ma.std(R),4)
figure()
hist(R.ravel(),bins=100,histtype='step',log=True)
show()

