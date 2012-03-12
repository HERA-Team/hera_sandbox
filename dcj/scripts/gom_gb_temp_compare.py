#!/usr/bin/env python
#
#  gom_gb_temp_compare.py
#  
#
#  Created by Danny Jacobs on 6/30/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
import numpy as n,pylab as p,datetime,gb_weather,ephem,traceback,sys
import aipy as a

file = open('concat_temp.txt','r'); data = [s.split('\t') for s in file.readlines()];file.close()
data.remove(['\n'])
T_GB=[]
T_CAB=[]
T_GOM = []
T_LOAD = []
T_RCVR = []

t = []
o = ephem.Observer()
t_last=0
for i,d in enumerate(data):
    if len(d)>9:
        mm,dd,yy = d[0].split('/')
        yy = '20'+yy
        hh,m,ss = d[1].split(':')
        yy,mm,dd,hh,m,ss = map(int,(yy,mm,dd,hh,m,ss))
        o.date = "%4d/%d/%d %d:%d:%d" %(yy,mm,dd,hh,m,ss)
        if o.date>t_last+30.0/(60*24):
            T_GOM.append(float(d[7])-273)
            T_RCVR.append(float(d[3])-273)
            T_LOAD.append(float(d[9])-273)    
            T_CAB.append(float(d[5])-273)
            try:
                t.append(datetime.datetime(yy,mm,dd,hh,m,ss))
            except:
                print traceback.print_exc()
                print o.date,i,d, 
                print  "%4d/%d/%d %d:%d:%d" %(yy,mm,dd,hh,m,ss)
                raise            
            T_GB.append(gb_weather.get_gb_temp(a.phs.ephem2juldate(o.date)+4.0/24))
            t_last = o.date
#        if i>100: break
#solve for linear relation between GB and GOM, etc
T_GOM,T_RCVR,T_LOAD,T_CAB,T_GB = map(n.array,(T_GOM,T_RCVR,T_LOAD,T_CAB,T_GB))
T_GB = T_GB.squeeze()
fig = p.figure()
print "offset, slope, rms"
for T in (T_GOM,T_RCVR,T_LOAD,T_CAB):
#    print T_GB.shape,T.shape
    D = n.vstack((n.ones_like(T_GB),T_GB)).transpose()
    G = n.linalg.lstsq(D,T)
#    print T.shape,D.shape,G#,G.shape,D.shape
    print G[0][0],G[0][1],n.sqrt(n.var(n.dot(D,G[0])-T))
    p.plot(t,n.dot(D,G[0]),'-',t,T,'.')
p.plot(t,T_GB,'-',linewidth=2)
p.xlabel('Local time (not UTC)')
p.ylabel('Temp [deg C]')
fig.autofmt_xdate()
p.show()
#for d in zip(T_GB,T_CAB,t):print d
#print i
sys.exit()
fig = p.figure()
p.plot(t,T_GB,'.',label="Green Bank standard")
p.plot(t,T_RCVR,'.',label="GoM Receiver")
p.plot(t,T_CAB,'.',label="GoM Cable")
p.plot(t,T_GOM,'.',label="GoM Balun")
p.plot(t,T_LOAD,'.',label="GoM Load")
p.xlabel('Local time (not UTC)')
p.ylabel('Temp [deg C]')
print min(t),max(t)
p.legend()
fig.autofmt_xdate()
p.show()