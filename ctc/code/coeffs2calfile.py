#!/usr/bin/env python

import numpy
import os, sys
import optparse

o = optparse.OptionParser()
o.add_option('--num', dest='num',type=int,default=3,
    help='number of coeffs##.npy files to use')
o.add_option('--refant1',dest='refant1',type=int,default=1,
    help='reference antenna 1 (not including main ref ant). Default is 1')
o.add_option('--refant2',dest='refant2',type=int,default=19,
    help='reference antenna 2 (not including main ref ant). Default is 19')
opts, args = o.parse_args(sys.argv[1:])

alldelays = numpy.zeros((7,16))
allgains = numpy.zeros((7,16))
refgains1 = []
refgains2 = []

for i in range(opts.num):
    file = 'coeffs'+str(("%02d" % int(i+1)))+'.npy'
    antpos,delays,gains = numpy.load(file)
    refantloc1 = numpy.where(antpos==opts.refant1)
    refantloc2 = numpy.where(antpos==opts.refant2)
    delays = numpy.load(file)[1]
    gains = numpy.load(file)[2]
    refgains1.append(gains[refantloc1][0])
    refgains2.append(gains[refantloc2][0])
    alldelays += delays
    allgains += gains
    #antpos[0] is first row
    
allgains -= opts.num-1
#print alldelays
#print allgains

gainref1 = numpy.mean(refgains1)
gainref2 = numpy.mean(refgains2) #averaging reference antenna gains

print 'Delays:'
for i in range(128):
    delay = alldelays[numpy.where(antpos==i)]
    if len(delay) != 0:
        print str(i)+': '+str(delay[0])+','
    else:
        print str(i)+': 0.0,'

print 'Gains:'
for i in range(128):
    gain = allgains[numpy.where(antpos==i)]
    if i == opts.refant1:
        print str(i) + ': '+str(gainref1)+ ','
    elif i == opts.refant2:
        print str(i) + ': '+str(gainref2)+','
    elif len(gain) != 0:
        print str(i)+': '+str(gain[0])+','
    else:
        print str(i)+': 0.0,'
