#!/usr/bin/env python
#
#  spectrum_fisher.py
#  
#
#  Created by Danny Jacobs on 4/17/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('--flim', default='0.13_0.18',
    help='Specify the spectrum limits. default=0.13_0.18')
o.add_option('--noise',default=1.0,
    help='Noise level in Jy [1]')
opts, args = o.parse_args(sys.argv[1:])

fmin,fmax = map(float,opts.flim.split('_'))
fmin,fmax = fmin/0.15,fmax/0.15
N = opts.noise


S = 10.0
B = -1.5

Gss = (fmax**(2*B+1) - fmin**(2*B+1))/(B)
GsB = S*(fmax**(2*B) - fmin**(2*B))
GBB = 2*(S*B)**2/(2*B-1)*(fmax**(2*B-1)-fmin**(2*B-1))

G = n.array([[Gss,GsB],[GsB,GBB]])/N
C = n.linalg.inv(G)
print C
V,W = n.linalg.eig(C)
print V
print W
