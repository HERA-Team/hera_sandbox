#! /usr/bin/env python
import sys, os, glob
try: f = sys.argv[1]
except(IndexError): sys.exit(0)
filelist = glob.glob('*.'+f.split('.')[-1])
filelist.sort()
try: i = filelist.index(f)
except(ValueError): sys.exit(0)
print ' '.join(filelist[i-1:i+2])


