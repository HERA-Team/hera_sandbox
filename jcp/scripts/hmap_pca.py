#!/usr/global/paper/bin/python

import numpy as n, pylab as p, aipy as a
import sys
#from PIL import Image

def convolve2d(a, b): return n.fft.ifft2(n.fft.fft2(a) * n.fft.fft2(b))

def get_data(filename):
    im = a.img.from_fits(filename)[0]
    return im

def make_scalar(d):
    #d = d.astype(n.float32)
    d = n.sqrt(n.sum(d**2, axis=2))
    d -= n.average(d)
    return d

maps = {}

for cnt, filename in enumerate(sys.argv[1:]):
    print 'Reading', filename
    d = get_data(filename)
    d = make_scalar(d)
    maps[cnt] = d

DIM = len(maps)
x = n.zeros((DIM,DIM), dtype=n.float)
for cnt1 in maps:
    b = maps[cnt1]
    for cnt2 in maps:
        c = n.average(b * maps[cnt2])
        print cnt1, cnt2, c
        x[cnt1,cnt2] = c

print x
u, s, v = n.linalg.svd(x)
u = u.transpose()
print u
for i in u:
    comp = 0
    for cnt, j in enumerate(i):
        comp += maps[cnt] * j
    p.imshow(comp)
    p.show()
        
#p.subplot(3,3,cnt+1)
#p.imshow(d)
#p.show()

#im = Image.fromarray(_dsum)
#im.save('test.jpg')


