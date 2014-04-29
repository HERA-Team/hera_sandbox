#! /usr/bin/env python

import sys, os, glob

try:
    f = sys.argv[1]
except(IndexError):
    sys.exit(0)

if f.endswith('/'):
    f = f[:-1]

path = os.path.dirname(f)

if len(path.strip()) != 0:
    filelist = sorted(glob.glob(path+'/*.'+f.split('.')[-1]))
else:
    filelist = sorted(glob.glob('*.'+f.split('.')[-1]))

try:
    i = filelist.index(f)
except(ValueError):
    sys.exit(0)

filelist = filelist[i-1:i+2]
print ' '.join(filelist)
