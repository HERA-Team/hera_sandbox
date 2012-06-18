#! /usr/bin/env python
import sys, os, glob
f = sys.argv[1]
filelist = glob.glob('*.'+f.split('.')[-1])
filelist.sort()
i = filelist.index(f)
print ' '.join(filelist[i-1:i+2])


