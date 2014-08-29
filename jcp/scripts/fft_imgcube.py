#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p, optparse, sys
import capo as C

o = optparse.OptionParser()
opts, args = o.parse_args(sys.argv[1:])

cube = []
for file in args:
    print "Reading", file
    map, kwds = a.img.from_fits(file)
    cube.append(map)
cube = n.array(cube)

#_cube = n.fft.ifft(cube,axis=0) #segfaults?!
_cube = n.fft.fft(cube,axis=0)/cube.shape[0]

delta = 4.8828125e-05 * 5 #5 channels per map
etas = n.fft.fftfreq(cube.shape[0],delta)
ks = C.pspec.dk_deta(C.pspec.f2z(.1525)) * etas
#print ks

print len(args), cube.shape, _cube.shape

for mode in xrange(cube.shape[0]):
    filename = 'psa747_null_k%.2f.bim.fits' % ks[mode]
    print 'Saving', filename
    #a.img.to_fits(filename,_cube[mode])
    a.img.to_fits(filename,n.abs(_cube[mode]),**kwds)

#p.imshow(n.abs(_cube[:,:,0]),interpolation='nearest')
#p.colorbar()
#p.show()
