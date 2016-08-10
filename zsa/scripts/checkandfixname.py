#! /usr/bin/env python

import sys,os,glob

try:f = sys.argv[1]
except(IndexError): sys.exit(0)

a,b,c,d = parts = f.split('.')
c = int(parts[2])
if f.endswith('/'): f = f[:-1]

if len(glob.glob('.'.join([a,b,str(c),d]))) != 0:
    fname = '.'.join([a,b,str(c),d])
if len(glob.glob('.'.join([a,b,str(c-1),d]))) != 0:
    fname = '.'.join([a,b,str(c-1),d])
if len(glob.glob('.'.join([a,b,str(c+1),d]))) != 0:
    fname = '.'.join([a,b,str(c+1),d])

print fname
#openfile = open(fixedfnames,'a')
#openfile.write(fbase+str(finalnum))
#openfile.close()

