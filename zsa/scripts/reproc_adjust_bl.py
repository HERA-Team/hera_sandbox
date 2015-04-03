#! /usr/bin/env python
import numpy as n 
import optparse, sys
import glob

o = optparse.OptionParser()

opts, args = o.parse_args(sys.argv[1:])
extension='*.uvcRREcAz'

def fsplit(name):
    name = name.split('.')
    t = name[1] + '.' + name[2]
    return t

frem = open('removed_from.txt','w')
for dir in args:
    print dir,
    files = glob.glob(dir + '/'+extension)
    if len(files) < 2: continue
    t0 = float(fsplit(files[0]))
    t1 = float(fsplit(files[1]))
    if t1-t0 < 0.04:
        print 'remove all files'
        frem.write(dir+'\n')
    else: print ' ' 
