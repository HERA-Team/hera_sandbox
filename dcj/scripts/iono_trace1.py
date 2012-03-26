#!/usr/bin/env python
#
#  iono_trace1.py
#  
#
#  Created by Danny Jacobs on 4/19/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

#load a fit trace for dra,ddec
#metrics
 # -make a sky image of, mag and gradient
 
def get_param_trace(root,ttype):
    """opens a a.fit.*param trace data set. consisting of root.pkl and root.txt
    then proceeds to mash it into a trace dict"""
#    prm_list = n.genfromtxt(root+'_prms.txt')
#    prm_keys = pickle.load(open(root+'_keys.pkl'))
    prm_trace = a.fit.reconstruct_prms(prm_list[0],prm_keys)
    prms = open(root+'_prms.txt')
    prm_keys =  open(root+'_keys.txt')
    #mash variables into arrays, assume its a 2 level param dict from a catalog.
    # (Why would I have different?)
    t = n.genfromtxt(root+'_times.txt',dtype={'names':['jd','lst'],'formats':[n.float,n.float]})
    times = None
    if root.endswith('fit'): 
        scores = n.loadtxt(root+'_scores.txt')
        prm_list = prm_list[scores[:,0]>0]
        t = t[scores[:,0]>0]
    for src,d in prm_trace.iteritems():
        for k in d.keys():
            if root.endswith('mdl'):
            #our bm/fit program currently misses the last time so chop those from the model too
#                print 'model hack',root
                prm_trace[src][k] = dict(zip(n.round(t[ttype],5),prm_list[:-1,prm_keys[src][k][0]]))
                if times is None: times = n.round(t[ttype][:-1],5)
            else:
                prm_trace[src][k] = dict(zip(n.round(t[ttype],5),prm_list[:,prm_keys[src][k][0]]))
                if times is None: times = n.round(t[ttype],5)
#                if k=='jys' and src=='J1143+221' and root.endswith('fit'):
#                     print prm_trace[src][k],len(n.round(t[ttype],5)),len(prm_list[:,prm_keys[src][k][0]])
#    if root.endswith('fit'):print prm_trace['J1143+221']['jys']
    return times,prm_trace
def interp_NaNs(Y):
    yleft = 0
    curnans = []
    for i,y in enumerate(Y):
        if n.isnan(y): curnans.append(i)
        else: 
            if len(curnans)>0 and not 0 in curnans:
                Y[curnans] = (yleft+y)/2.
            elif len(curnans)>0 and 0 in curnans:
                Y[curnans] = y
            curnans = []
            yleft = y
        if n.isnan(y) and i+1==len(Y):
            Y[curnans] = yleft
    return Y
    
    
    
root = args[0]
fit_times,fit_params = get_param_trace(root+'/fit','lst')    
fit_scores = loadtxt(root+'/fit_scores.txt')

