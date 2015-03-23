#!/usr/global/paper/bin/python
import aipy as a, numpy as n, optparse, sys, pylab as p 

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

freqs, pixvals = [],[]
for file in args:
    chan = (float(file.split('_')[2][1:]) + float(file.split('_')[3]))/2.
    #chan = (float(file.split('_')[1][1:]) + float(file.split('_')[2]))/2.
    freq = chan * (.1/2048) + .1
    d,kwds = a.img.from_fits(file)
    d = d[::-1,]
    #xpx,ypx = 0.9725,0.7825 #coords from plot_img
    xpx,ypx = 576,421
    #xpx,ypx = 576+1,421+1
    xpx,ypx = 586,511
   
    freqs.append(freq) 
    pixvals.append(d[xpx,ypx])

inds = n.argsort(freqs)
freqs, pixvals = n.take(freqs,inds), n.take(pixvals,inds)

#print n.mean(pixvals),n.mean(pixvals[0:len(pixvals)/2])
#print n.std(pixvals),n.std(pixvals[0:len(pixvals)/2])
print n.sum(pixvals)

p.plot(freqs,pixvals)
p.show()
