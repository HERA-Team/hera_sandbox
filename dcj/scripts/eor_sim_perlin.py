#! /usr/bin/env python
from pylab import *
import noise

RAsize = 4096*2
Decsize = 1024
T = n.zeros((RAsize,Decsize))
subplot(211)
for i in xrange(RAsize):
    for j in xrange(Decsize):
        T[i,j] = noise.pnoise2(i/16./octaves,j/16./octaves,octaves,lacunarity=2)
imshow(T.T)


for i in xrange(RAsize):
    for j in xrange(Decsize):
        T[i,j] = noise.pnoise2(i/16./octaves,j/16./octaves,octaves,lacunarity=2,base=16)
T[T<0] = 0.
subplot(212)
imshow(T.T)
savefig('eor_sim_perlin.png')
