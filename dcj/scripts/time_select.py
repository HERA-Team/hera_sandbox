#!/usr/bin/env python
#
#  time_select.py
#  
#
#  Created by Danny Jacobs on 5/20/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,os

o = optparse.OptionParser()
o.set_usage('time_select.py -t <start>_<stop> -u lst <files>')
a.scripting.add_standard_options(o, cal=True)
o.add_option('-t', 
    help='Time range (<t1>_<t2>) in units specified by -u. Use : to indicate either end.')
o.add_option('-u',default='jd',
    help='Time units. (jd,lst[decimal hrs]) default jd')
opts, args = o.parse_args(sys.argv[1:])
aa = a.cal.get_aa(opts.cal,0.1,0.0001,1)
curtime = 0
if not opts.u in ['jd','lst']:
    raise(ValueError, "time unit %s not supported"%opts.u)
try:
    t1,t2 = opts.t.split('_')
    if t1==':' and opts.u=='jd': t1=-n.Inf
    elif t1==':' and opts.u=='lst':t1=0
    else: t1 = float(t1)
    if t2==':' and opts.u=='jd': t2=n.Inf
    if t2==':' and opts.u=='lst': t2=2*n.pi
    else: t2 = float(t2)
    if opts.u=='lst': t1,t2 = t1*n.pi/12,t2*n.pi/12
except(ValueError):
     print "Error parsing %s. Please enter a start and a stop time"%opts.t
if opts.u=='lst' and (t1<0 or t1 > 2*n.pi): 
    raise(ValueError, "Values %s out of range for lst"%opts.t)

for filename in args:
    print filename, '->', filename+'t'
    if os.path.exists(filename+'t'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'t', status='new') 
    
    def mfunc(uv,p,d):
        crd,t,(i,j) = p
        if opts.u=='lst':t=uv['lst']
        if t>=t1 and t<t2: return p,d
        else: return p,None
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,
        append2hist='TIME_SELECT:   time range: %s,  units:%s\n'%(opts.t,opts.u))
    del(uvo)