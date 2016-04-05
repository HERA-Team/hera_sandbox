#! /usr/bin/env python
import numpy as n, pylab as p
import sys, aipy as a,re
i,j = 64,49
mybl=str(a.miriad.ij2bl(i,j))
jdjump = 2456949.
times = []
print mybl
def file2jd(zenuv):
    return float(re.findall(r'\d+\.\d+', zenuv)[0])
for filename in sys.argv[1:]:
    print 'Reading', filename
    try:
        npz = n.load(filename)
    except:
        print "    failed to load"
    if file2jd(filename)>jdjump:c='orange'
    else: c='k'
    for k in npz.files:
	if not k.endswith(mybl): continue
        if k.startswith('t'): continue
        p.figure(1)
	p.subplot(311)
	p.plot(npz['t'+k]*12/n.pi,n.abs(npz[k]),color=c)
        p.subplot(312)
        p.plot(npz['t'+k]*12/n.pi,n.angle(npz[k]),color=c)
        p.subplot(313)
        p.plot(npz['t'+k]*12/n.pi, n.real(npz[k]),color=c)
p.figure(1)
p.subplot(311)
p.title('color changes at {jdjump}, baseline:{i}_{j}'.format(jdjump=jdjump,i=i,j=j))

p.show()
