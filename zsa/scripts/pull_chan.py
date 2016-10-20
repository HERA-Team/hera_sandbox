#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import glob, sys, os,optparse

o = optparse.OptionParser()
o.add_option('-c', "--chan", help="Channels to pull",type=str)
o.add_option('--bls',help='Baselines, comma delim eg "41_64,65_66" default is all bls of type 0N,2E in the psa128 array', default='auto')
o.add_option('-p','--pol',help='polarization. For now, you must pick one of: xx,yy,xy or yx')
opts, args = o.parse_args()
filelist = sys.argv[1:]
bls=opts.bls

print "using bls:",bls
for f in args:
    npz = {}
    print 'Reading', f
    outfile = os.path.basename(f) + '.npz'
    if os.path.exists(outfile):
        print outfile, "exists, skipping"
    uv = a.miriad.UV(f)
    CHAN = a.scripting.parse_chans(opts.chan, uv['nchan'])
    a.scripting.uv_selector(uv, ants=bls, pol_str=opts.pol)
    for (crd,t,(i,j)),d,_f in uv.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        sbl = str(bl)
        if len(d)<n.max(CHAN): print "ERROR Some Channels not found in data range: 0",len(d)
        npz[sbl] = npz.get(sbl,[]) + [d[CHAN]]
        npz['t'+sbl] = npz.get('t'+sbl,[]) + [uv['lst']]
    print 'Writing', outfile
    n.savez(outfile, **npz)
        

#dsort = {}
#for k in data:
#    print k, len(data[k])
#    i = n.argsort(time[k])
#    dsort[str(k)] = n.array(data[k])[i]
#    dsort['t'+str(k)] = n.array(time[k])[i]
    


