#!/usr/global/paper/bin/python

import aipy as a, ephem as e, numpy as np, sys, optparse, os, math as m


o = optparse.OptionParser()
o.set_usage('testbed.py [options]')
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, src=True, prms=True)

opts, args = o.parse_args(sys.argv[1:])
opts.ant += ',cross'

uv = a.miriad.UV('/data1/paper/arp/pgb015/zen.2455015.72815.uvc')
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])

prms, prm_dict, shkeys = {}, {}, []
pd = a.scripting.parse_prms(opts.prms)
for k in pd:
    if prms.has_key(k): print 'yippeeeeeeeeeeeeeeeeeeeeeeeeeeeee'
    else: prms[k] = pd[k]
for prm in prms: prm_dict[prm] = prms[prm].keys()
start_prms = aa.get_params(prm_dict)
start_prms.update(cat.get_params(prm_dict))

for obj in start_prms:
    for prm in start_prms[obj]:
        if prms[obj][prm][0] !=None:
            start_prms[obj][prm] = prms[obj][prm][0]

prm_list, key_list = a.fit.flatten_prms(start_prms)

fit_fit = None
mfq = cat.get('mfreq')
dbuf = None

print prms

def fit_func(prms, filelist, decimate, decphs):
    print 'you are now in functionland'
