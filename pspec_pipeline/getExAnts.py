#! /usr/bin/env python

import sys, optparse, re
o = optparse.OptionParser()
o.add_option('-f','--badAntsFile',dest='badantfile',help='path to text file organized as: JD<space>comma,separated,bad,antnums<CR>')
opts,args = o.parse_args(sys.argv[1:])

assert(len(args)==1) #one file at a time only (for now)

def file2jd(zenuv):
    return float(re.findall(r'\d+\.\d+', zenuv)[0])

jd = int(file2jd(args[0]))

#this could probably be done more efficiently
with open(opts.badantfile) as search:
    for line in search:
        if line.startswith(str(jd)):
            print line.split(' ')[1]
