#! /usr/bin/env python
import numpy as np, sys, matplotlib.pyplot as plt, optparse
o = optparse.OptionParser()
o.set_usage('read_flags.py [options] *xx.txt')
o.set_description('parse text output of meanVij.py')
o.add_option('--epochJD', default=None, help='JD _after which_ there is an epoch boundary.') #XXX this only works if there are two Epochs
o.add_option('-v','--verbose',action='store_true',help='toggle verbosity')

opts, files = o.parse_args(sys.argv[1:])

if not opts.epochJD is None: 
    opts.epochJD = int(opts.epochJD)
    go = True
    
grid = np.zeros((len(files),112),dtype='int')
JDs = []
start,end = 0,len(files)

for i,f in enumerate(files):
    with open(f) as f: lines = f.readlines()
    jd,bas = lines[0].rstrip().split('    ')
    jd,bas = int(jd), map(int,bas.split(','))
    JDs.append(jd)
    
    #find epoch bound in a silly way. must be an array manip way to do this.
    if not opts.epochJD is None and go:
        if jd>=opts.epochJD:
            end = i
            go = False
    
    if opts.verbose: print jd,bas
    for ba in bas: grid[i,ba] = 1

JDs=np.array(JDs)    

if not opts.epochJD is None:
    c1 = np.sum(grid[start:end,:],axis=0)
    c2 = np.sum(grid[end:,:],axis=0)
    print 'EPOCH 1'
    print np.argsort(c1)[::-1]
    print '======='
    print 'EPOCH 2'
    print np.argsort(c2)[::-1]

else:
    c = np.sum(grid,axis=0)
    print np.argsort(c)[::-1]
