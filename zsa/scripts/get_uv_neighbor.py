#! /usr/bin/env python
import sys, os, glob
from numpy import abs,diff
# Make sure we have an argument
try: f = sys.argv[1]
except(IndexError): sys.exit(0)
if f.endswith('/'): f = f[:-1] #the script fucks up if you give it a directory...
path = os.path.dirname(f)
if len(path.strip()) != 0: filelist = glob.glob(path+'/*.'+f.split('.')[-1])
else: filelist = glob.glob('*.'+f.split('.')[-1])
filelist.sort()
try: i = filelist.index(f)
except(ValueError): sys.exit(0)

#Check that files are actually sequential (dfm):
filelist = filelist[i-1:i+2]
filedates= [float(f.split('.')[1]+'.'+f.split('.')[2]) for f in filelist]
if (abs(diff(filedates)) > 20./(24.*60.)).any():
    sys.exit(0) #if files are more than 20 minutes apart, break.
else: print ' '.join(filelist)
