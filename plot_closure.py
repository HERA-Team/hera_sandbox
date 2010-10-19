#!/usr/bin/env python
#
#  plot_closure.py
#  
#
#  Created by Danny Jacobs on 6/2/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,time,resource
from pylab import *
import logging
from matplotlib.widgets import RadioButtons
#from dcj_plot import plot_traces
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('plot_closure')

I = n.complex(0,1)

o = optparse.OptionParser()
a.scripting.add_standard_options(o, pol=True)
o.set_usage('use the antenna option to specifiy a single triangle over which to compute phase closure')
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('-a',help='specify two antennae between which to compute phase closure [eg 1,2]')
o.add_option('--via',help='Specify antennae number range. eg 0_31')
o.add_option('-c','--chan',help='Supply and optional single channel range. eg -c100_300')
opts, args = o.parse_args(sys.argv[1:])
if not opts.chan is None:
    print opts.chan
    chans  = map(int,opts.chan.split('_'))
    chans = n.arange(chans[0],chans[1])
else:
    chans = None
short = '_'.join(opts.a.split(','))
ants =  map(int,opts.a.split(','))
vias = map(int,opts.via.split(','))
vias = range(vias[0],vias[-1]+1)
#vias = list(set(vias).difference(set(ants)))
longs = {}
for via in vias:
    longs[via]=['%s_%i'%(ants[0],via),'%i_%s'%(via,ants[-1])]
print longs
phases = {}
ttype = n.dtype([('jd',float),('lst',float),('index',int)])
ind = 0
times = []
bls = []
for filename in args:
    uv = a.miriad.UV(filename)
    print "reading: ",filename

    a.scripting.uv_selector(uv, short+','+','.join([','.join(c) for c in longs.values()])+',cross', opts.pol)
#    print short+','+','.join([','.join(c) for c in longs.values()])+',cross'

    curtime = 0
    for p,d in uv.all():
        uvw,t,(i,j) = p

        if not chans is None:
            d = d.take(chans)
        bl = '%s_%s'%(i,j)
        if t != curtime:
            phases[t] = {}
            curtime = t
            times.append((t,uv['lst']*12/n.pi,ind))
            ind += 1
        for ant,long in longs.iteritems():
            if not phases[t].has_key(ant):
                phases[t][ant] = []
            if bl in long:
                phases[t][ant].append(-n.angle(d))
                if ind==1: bls.append(bl)
#                print '+',
            elif '_'.join(bl.split('_')[::-1]) in long:
                phases[t][ant].append(-n.angle(d))
                if ind==1: bls.append(bl)
#                print '-',
        if bl == short:
            try: phases[t][ant].append(-n.angle(d))
            except(KeyError):phases[t][ant] = [-n.angle(d)]
            if ind==1: bls.append(bl)
#        elif bl[::-1] == short:
#            try: phases[t][ant].append(n.angle(d))
#            except(KeyError):phases[t][ant] = [n.angle(d)]
print resource.getrusage(resource.RUSAGE_SELF)[8]

times = n.array(times,dtype=ttype)
#print bls
#figure()
#for i,y in enumerate(phases[times['jd'][0]][8]):
#    plot(y,label=bls[i])
#plot(range(len(y)),n.sum(phases[times['jd'][0]][8],axis=0),'--',label='sum')
#legend()
#show()
#sys.exit()
closures = {} 
for t,ps in phases.iteritems():
    for ant,P in ps.iteritems():
        if not closures.has_key(ant):
            closures[ant] = {}
        if ant in ants:
            closures[ant][t] = n.zeros_like(d)
        else:
            closures[ant][t] = n.sum(n.array(P),axis=0)
            

jd =sort(times['jd'])
ants = list(sort(closures.keys()))
if len(ants)==1:
    D = n.array([closures[ants[0]][t] for t in jd])
    figure()
    C = 0
    subplot(211)
    for i in range(3):
        plot([phases[t][ants[0]][i][C] for t in jd],label=bls[i])
    plot(D[:,C],'k',label='sum')
    legend()
    subplot(212)
    for i in range(3):
        d = n.exp(-I*n.array(phases[jd[0]][ants[0]][i]))
        id = n.fft.ifft(d)
        id = n.concatenate([id[id.shape[0]/2:], id[:id.shape[0]/2]], 
                axis=0)
        plot(id,label=bls[i])
    id = n.fft.ifft(n.exp(-I*D[0,:]))
    id = n.concatenate([id[id.shape[0]/2:], id[:id.shape[0]/2]], 
            axis=0)    
    plot(id,'k',label='sum')
    figure()

    subplot(311)
    imshow(D,aspect='auto',extent=(chans[0],chans[-1],D.shape[0],0))
    colorbar(orientation='horizontal',shrink=0.75,ticks=(n.ceil(n.min(D)),0,n.floor(n.max(D))))
    subplot(312)
    hist(n.average(D,axis=1),bins=n.sqrt(D.shape[1]))
    xlabel('phase [r]')
    subplot(313)
    hist(n.average(D,axis=0),bins=n.sqrt(D.shape[0])/10)
    show()
if len(ants)>1:
    fclosure = n.real(n.array([n.average([closures[ant][t] for t in jd],axis=0) for ant in ants]))
    fclosure_std = [n.std([closures[ant][t] for t in jd],axis=0) for ant in ants]
    tclosure = n.array([n.average([closures[ant][t] for t in jd],axis=1) for ant in ants])
    tclosure = n.real(tclosure.transpose())
    tclosure_std = n.array([n.std([closures[ant][t] for t in jd],axis=1) for ant in ants]).transpose()


    #sys.exit()
    fig0 = figure()
    fig0.canvas.manager.set_window_title('plot_closure waterfalls std '+time.strftime("%Y-%m-%d %H:%M:%S"))

    subplot(121)
    twf = imshow(tclosure_std,aspect='auto',interpolation='nearest',extent=(n.min(vias),n.max(vias),n.max(times['index']),n.min(times['index'])))
    twfa = twf.get_axes()
    colorbar(orientation='horizontal',shrink=0.75,ticks=(0,n.int(n.max(tclosure))))
    print 'max tclosure = ',round(n.max(tclosure),1)
    title('time-baseline closure')
    xlabel('baseline from %s to %s via ..'%(opts.a.split(',')[0],opts.a.split(',')[-1]))
    ylabel('index')

    subplot(122)
    fwf = imshow(fclosure_std,aspect='auto',interpolation='nearest',extent=(0,d.shape[0],n.max(vias),n.min(vias)))
    ylabel('baseline from %s to %s via ..'%(opts.a.split(',')[0],opts.a.split(',')[-1]))
    xlabel('channel')
    title('frequency-baseline closure')
    print 'max fclosure = ',round(n.max(fclosure),1)
    colorbar(orientation='horizontal',shrink=0.75,ticks=(0,n.int(n.max(fclosure))))

    subplots_adjust(bottom=0.3)
    suptitle(opts.a)    
    
    fig1 = figure()
    fig1.canvas.manager.set_window_title('plot_closure waterfalls '+time.strftime("%Y-%m-%d %H:%M:%S"))

    subplot(121)

    twf = imshow(tclosure,aspect='auto',interpolation='nearest',extent=(n.min(vias),n.max(vias),n.max(times['index']),n.min(times['index'])))
    twfa = twf.get_axes()
    colorbar(orientation='horizontal',shrink=0.75,ticks=(0,n.floor(n.max(tclosure))))

    title('time-baseline closure')
    xlabel('baseline from %s to %s via ..'%(opts.a.split(',')[0],opts.a.split(',')[-1]))
    ylabel('index')

    subplot(122)
    fwf = imshow(fclosure,aspect='auto',interpolation='nearest',extent=(0,d.shape[0],n.max(vias),n.min(vias)))
    ylabel('baseline from %s to %s via ..'%(opts.a.split(',')[0],opts.a.split(',')[-1]))
    xlabel('channel')
    title('frequency-baseline closure')

    colorbar(orientation='horizontal',shrink=0.75,ticks=(0,n.floor(n.max(fclosure))))

    subplots_adjust(bottom=0.3)
    suptitle(opts.a)
    bax = axes([.3,0.025,0.15,0.15],axisbg='lightgoldenrodyellow')
    rtime = RadioButtons(bax,('jd','lst','index'))
    def tfunc(label):
        tstart,tstop = n.min(times[label]),n.max(times[label])
        twf.set_extent((n.min(vias),n.max(vias),tstop,tstart))
        twfa.set_ylabel(label)
        #fwf.set_extent((n.min(vias),n.max(vias),tstop,tstart))
        draw()
    rtime.on_clicked(tfunc)    
    draw()
    fig2 = figure()
    fig2.canvas.manager.set_window_title('plot_closure averages '+time.strftime("%Y-%m-%d %H:%M:%S"))
    subplot(211)
    plot(n.average(tclosure,axis=1))
    xlabel('time')
    subplot(212)
    plot(n.average(fclosure,axis=0))
    xlabel('frequency')
    show()

        