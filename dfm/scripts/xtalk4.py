#!/usr/bin/env python
"""
Remove cross-talk by subtracting long-time (nightly) average from the data.
"""
import aipy as a
import numpy as np
import optparse,sys,os
from glob import glob

o = optparse.OptionParser()
o.set_usage('xtalk4.py *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

#args = glob('array0/zen.*.uvA0rcb')
#args.sort()

xtalk,cnt = {},{}

for uvfile in args:
    print '='*(len(uvfile)+8)
    print 'Reading %s'%uvfile
    uv = a.miriad.UV(uvfile)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        pol = a.miriad.pol2str[uv['pol']]
        bl = (i,j)
        if i == j: continue
        if not pol in xtalk.keys(): 
            xtalk[pol] = {}
            cnt[pol] = {}
        if not bl in xtalk[pol].keys():
            xtalk[pol][bl] = np.zeros(uv['nchan'],dtype=complex)
            cnt[pol][bl] = np.zeros(uv['nchan'])
        xtalk[pol][bl] += np.where(f,0,d)
        cnt[pol][bl] += np.where(f,0,1)
    del uv
print '='*(len(args[-1])+8)

for pol in xtalk:
    for bl in xtalk[pol]:
        xtalk[pol][bl] /= cnt[pol][bl]

for infile in args:
    outfile = infile+'X'
    print infile,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping...'
        continue
    uvi = a.miriad.UV(infile)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)

    def mfunc(uv,p,d):
        uvw,t,(i,j) = p
        pol = a.miriad.pol2str[uv['pol']]
        bl = (i,j)
        if i != j: d -= xtalk[pol][bl]
        return p,d

    uvo.pipe(uvi,mfunc=mfunc,append2hist='XTALK4 \n')
