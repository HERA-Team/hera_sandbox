__author__ = 'yunfanzhang'

import numpy as n
import random

def cross_mult(data1, data2):
    return n.multiply(n.conjugate(data1), data2)

#Bootstrap over a 2d array
def bootstrap(B, data):
    boot = []
    for b in range(B):
        temps = n.array([],dtype='complex64')
        for i in range(len(data[0])):    #for each frequency channel
            temps = n.append(temps,0.+0.j)
            for j in range(len(data)):      #bootstrap over all the time samples
                choi = random.choice(data)   #choi contains all frequency channels of a time sample
                temps[i] = temps[i] + choi[i]
        boot.append(temps/float(len(data)))
    boot = n.array(boot).transpose()
    print "shape of bootstrap samples", boot.shape
    bootmean,booterr = [],[]
    for ch in range(len(boot)):
        mean = n.sum(boot[ch])/float(B)
        sig = n.sqrt(n.sum((boot[ch]-mean)**2)/(B-1.))
        bootmean.append(mean)
        booterr.append(sig)

    return bootmean, booterr
#p.hist(boot[5])
#p.show()

#simple test
#x = n.array([[1, 2, 3], [4, 5, 6],[7,8,9],[10,11,12],[13,14,15]])
#bm, be = boot_simple.bootstrap(100, x)