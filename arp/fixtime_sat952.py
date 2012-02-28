#! /usr/bin/env python

import aipy, sys, os, numpy, optparse

o = optparse.OptionParser()
o.set_usage('fix_time.py [options] *.uv')
o.set_description(__doc__)

opts, args = o.parse_args(sys.argv[1:])

last_time = 0
index_cnt = 1
bl_cnt = 0

def mfunc(uv, p, d, f):
    uvw, t, (i,j) = p
    global last_time, index_cnt, bl_cnt
    if not (last_time == t):
        #print t, last_time, index_cnt
        index_cnt = 1
        bl_cnt = 0
        last_time = t
    else:
        last_time = t
    	bl_cnt += 1
    	if bl_cnt == 10:
    	    index_cnt += 1
    	    bl_cnt = 0
        t = t + 0.000001 * index_cnt
        #print "%f"%(t)
    p = (uvw,t,(i,j))
    return p, d, f

for filename in sys.argv[1:]:
    print filename
    if os.path.exists(filename+'.t'):
        print '    File exists... skipping.'
        continue
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(filename+'.t', status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='Fix time resolution\n')
    del(uvo)