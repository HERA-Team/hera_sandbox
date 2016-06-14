#! /usr/bin/env python
"""
Plot antenna locations in a three column file
usage
plot_MWA_ants.py antfile.txt

cat antfile.txt
X Y Z
X Y Z

"""
from pylab import *
import numpy as n,sys
def i2a(i):
    #convert from MWA to unit ant index
    tens = int(i/8)+1
    ones = i%8+1
    return tens*10+ones
def a2i(a):
    #convert from unit to MWA ant index
    eights = int(a/10)-1
    ones = a%10
    return eights*8+ones
A = n.loadtxt(sys.argv[1])
Z = n.mean(A,axis=0)
#Z = A[:,0]
Z.shape = (1,) + Z.shape
A -= Z
labels = [str(i2a(i)) for i in range(0,128)]
#for i in xrange(128):
#    print i,i2a(i)
scatter(A[:,0],A[:,1],marker='x')
for i,l in enumerate(labels):
    text(A[i,0],A[i,1],l)
    print i,l,A[i,0],A[i,1]
xlim([-1600,1600])
ylim([-1600,1600])
show()

