#! /usr/bin/env python
import numpy as n, pylab as p, optparse, sys, ast, aipy as a

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

file = args[0]
data = open(file,'r').read().split('\n')

times = []
cas = {'flux':[],'dra':[],'ddec':[]}
cyg = {'flux':[],'dra':[],'ddec':[]}

data = data[:-1]

for ind,int in enumerate(data):
    if ind % 3 == 0:
        times.append(float(int))
    if ind % 3 == 1:
        continue
    if ind % 3 == 2:
        dict = ast.literal_eval(int)
        cas['flux'].append(dict['cas']['jys'])
        cas['dra'].append(dict['cas']['dra']/a.const.arcmin)
        cas['ddec'].append(dict['cas']['ddec']/a.const.arcmin)
        cyg['flux'].append(dict['cyg']['jys'])
        cyg['dra'].append(dict['cyg']['dra']/a.const.arcmin)
        cyg['ddec'].append(dict['cyg']['ddec']/a.const.arcmin)

#times = times[:-1]

p.subplot(231)
p.title('cas')
p.plot(times,cas['flux'])
p.ylabel('jys')
p.subplot(232)
p.plot(times,cas['dra'])
p.ylabel('dra')
p.subplot(233)
p.plot(times,cas['ddec'])
p.ylabel('ddec')
p.subplot(234)
p.title('cyg')
p.plot(times,cyg['flux'])
p.ylabel('jys')
p.subplot(235)
p.plot(times,cyg['dra'])
p.ylabel('dra')
p.subplot(236)
p.plot(times,cyg['ddec'])
p.ylabel('ddec')
p.show()
