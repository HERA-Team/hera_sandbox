#! /usr/bin/env python

import aipy as a, sys, optparse, os

rewire_run1 = { 
    0: 0,
    1: 1,
    2: 2,
    3: 3,
    4: 4,
    5: 5,
    6: 6,
    7: 8,
    8: 9,
    9: 10,
    10: 11,
    11: 12,
    12: 13,
    13: 14,
    14: 15,
    15: 16,
    16: 17,
    17: 18,
    18: 19,
    19: 20,
    20: 21,
    21: 22,
    22: 23,
    23: 24,
    24: 25,
    25: 26,
    26: 27,
    27: 28,
    28: 29,
    29: 30,
    30: 31,
    31: 32
}

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
opts,args = o.parse_args(sys.argv[1:])

def mfunc(uv, p, d):
    crd, t, (i, j) = p
    i,j = rewire_run1[i], rewire_run1[j]
    p = crd, t, (i, j)
    return p, d
 
for filename in args:
    print filename, '->', filename+'C'
    uvi = a.miriad.UV(filename)
    if os.path.exists(filename+'C'):
        print '    File exists... skipping.'
        continue
    uvo = a.miriad.UV(filename+'C',status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,append2hist='Corrected antenna assignment.\n')
