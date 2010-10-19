#!/usr/bin/env python
#
#  plot_closure_all.py
#  
#
#  Created by Danny Jacobs on 6/6/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,time
from pylab import *

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,pol=True)
o.add_option('-a','--refant',type='int',help='Reference antenna')
o.add_option('--via',help='Specify antenna range for vias and directs.')
o.add_option('-c','--chan',help='An option single channel range')
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])



def word_wrap(string, width=80, ind1=0, ind2=0, prefix=''):
    """ word wrapping function.
        string: the string to wrap
        width: the column number to wrap at
        prefix: prefix each line with this string (goes before any indentation)
        ind1: number of characters to indent the first line
        ind2: number of characters to indent the rest of the lines
    """
    string = prefix + ind1*" " + string
    newstring = ""
    if len(string) > width:
        while True:
            # find position of nearest whitespace char to the left of "width"
            marker = width-1
            while not string[marker].isspace():
                marker = marker - 1

            # remove line from original string and add it to the new string
            newline = string[0:marker] + "\n"
            newstring = newstring + newline
            string = prefix + ind2*" " + string[marker+1:]

            # break out of loop when finished
            if len(string) <= width:
                break
    
    return newstring + string
    
    
if not opts.chan is None:
    print opts.chan
    chans  = map(int,opts.chan.split('_'))
    chans = n.arange(chans[0],chans[1])
else:
    chans = None
A = opts.refant    
Bs = map(int,opts.via.split(','))
Bs = n.arange(Bs[0],Bs[-1]+1).astype(n.int)
Cs = Bs.copy()


def compute_closure(D,A,Bs,Cs):
    """
    Given a dictionary of baseline:spectrum pairs, compute average phase 
    closure between A and all of Cs averaged over all of Bs.
    """
    C_closure={}
    for C in Cs:
        B_closure = {}
        for B in Bs:
            if A==B or B==C or A==C:continue
            try:B_closure[B] = n.angle(D['%i_%i'%(A,B)])
            except(KeyError):B_closure[B] = - n.angle(D['%i_%i'%(B,A)])
            try: B_closure[B] += n.angle(D['%i_%i'%(B,C)])
            except(KeyError): B_closure[B] += -n.angle(D['%i_%i'%(C,B)])
            try: B_closure[B] += -n.angle(D['%i_%i'%(A,C)])
            except(KeyError): B_closure[B] += n.angle(D['%i_%i'%(C,A)])
        if A!=C:C_closure[C] = n.average(B_closure.values(),axis=0)
        else: C_closure[C] = n.zeros_like(D[D.keys()[0]])
    return C_closure

times =[]
ind =0
closures = {}
for filename in args:
    uv = a.miriad.UV(filename)
    print "reading: ",filename
    a.scripting.uv_selector(uv,'cross')
    if not opts.cal is None:
        aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['freq'],uv['nchan'])
    
    curtime = 0
    for p,d in uv.all():
        uvw,t,(i,j) = p
        if t-curtime>2e5:
            D = {}
            told = t
            curtime = t
        if not chans is None:
            d = d.take(chans)
        bl = '%i_%i'%(i,j)
        D[bl] = d
        if t!= curtime and t-curtime<2e5:
            if not opts.cal is None:
                aa.set_jultime(t)
            closures[curtime] = compute_closure(D,A,Bs,Cs)
            del(D)
            D = {}
            times.append((curtime,uv['lst']*12/n.pi,ind))
            ind+=1
            curtime = t
        
ttype = n.dtype([('jd',float),('lst',float),('index',int)])
times = n.array(times,dtype=ttype)
jd = n.sort(times['jd'])
Cs = n.sort(Cs)

#print len(jd), len(closures.keys()),n.max(jd),n.max(closures.keys())
#print [n.average([closures[t][ant] for t in jd],axis=0) for ant in Cs]
fclosure = n.real(n.array([n.average([closures[t][ant] for t in jd],axis=0) for ant in Cs]))
fclosure_std = [n.std([closures[t][ant] for t in jd],axis=0) for ant in Cs]
tclosure = n.array([n.average([closures[t][ant] for t in jd],axis=1) for ant in Cs])
tclosure = n.real(tclosure.transpose())
tclosure_std = n.array([n.std([closures[t][ant] for t in jd],axis=1) for ant in Cs]).transpose()    



fig0 = figure()

fig0.canvas.manager.set_window_title('plot_closure waterfalls '+time.strftime("%Y-%m-%d %H:%M:%S"))

subplot(121)
twf = imshow(tclosure,aspect='auto',interpolation='nearest',extent=(n.min(Cs),n.max(Cs)+1,n.max(times['index']),n.min(times['index'])))
twfa = twf.get_axes()
colorbar(orientation='horizontal',shrink=0.75,ticks=(n.int(n.min(tclosure)),0,n.int(n.max(tclosure))))
print 'max tclosure = ',round(n.max(tclosure),1)
title('time-baseline closure')
xlabel('average closure from %s to ...'%(opts.refant))
ylabel('index')

subplot(122)
fwf = imshow(fclosure,aspect='auto',interpolation='nearest',extent=(0,d.shape[0],n.max(Cs)+1,n.min(Cs)))
ylabel('average closure from %s to ...'%(opts.refant))
xlabel('channel')
title('frequency-baseline closure')
print 'max fclosure = ',round(n.max(fclosure),1)
colorbar(orientation='horizontal',shrink=0.75,ticks=(n.int(n.min(tclosure)),0,n.int(n.max(fclosure))))

subplots_adjust(bottom=0.3)
suptitle(opts.refant)    
bax = axes([.01,0.01,0.8,0.15],axisbg='lightgoldenrodyellow')
bax.set_axis_off()
text(0,0,word_wrap(' '.join(sys.argv),width=100))
#text(0, 0, r"\parbox[b]{4 in}{really really really really really really long line}", va="top")
show()