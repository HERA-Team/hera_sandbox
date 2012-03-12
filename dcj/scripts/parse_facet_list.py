#!/usr/bin/env python
#
#  parse_facet_list.py
#  
#
#  Created by Danny Jacobs on 9/4/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,re

o = optparse.OptionParser()
opts, args = o.parse_args(sys.argv[1:])

#M = re.compile('.*\(([-+]?[0-9]*\.?[0-9]+), ([-+]?[0-9]*\.?[0-9]+) in deg\).*')
M = re.compile('RA=(.*)\sDEC=(.*)$')
filename = args[0]
File = open(filename,'r')
srcs = []
for line in File.readlines():
    M.search(File.read())
    if not M is None:
        srcs.append(M.groups()[0]+M.groups()[1])
print ' '.join(srcs)
