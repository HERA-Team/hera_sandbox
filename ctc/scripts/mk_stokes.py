#! /usr/bin/env python

import aipy as a, numpy as n
import capo as C
import optparse, sys, os
import glob

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
o.add_option('--stokes', help='Comma-delimited list of Stokes parameters to form (e.g. I,Q,U,V).')
opts,args = o.parse_args(sys.argv[1:])


POL_WGTS = {
    'I': {'xx': 1. , 'yy': 1. },
    'Q': {'xx': 1. , 'yy':-1. },
    'U': {'xy': 1. , 'yx': 1. },
    'V': {'xy':-1.j, 'yx': 1.j},
}

def mfunc(uv,p,d,f):
    if p[2] in d2.keys():
        data2 = d2[p[2]][POL_WGTS[s].keys()[1]]
        flags2 = f2[p[2]][POL_WGTS[s].keys()[1]]
        data2 = data2[n.argmin(n.abs(t2['times']-p[1]))] #time index to grab from d2
        flags2 = flags2[n.argmin(n.abs(t2['times']-p[1]))]
        d = ((d*POL_WGTS[s][POL_WGTS[s].keys()[0]]) + (data2*POL_WGTS[s][POL_WGTS[s].keys()[1]]))/2.0
        f = f + flags2
        d = n.ma.masked_array(d,f)
        return p,d,f
    else:
        return p,None,None

stokes = opts.stokes.split(',')

jds = [] #unique files
for f in args:
    jds.append(f.split('.')[1]+'.'+f.split('.')[2])
    file_ending = f.split('.')[-1]
jds = sorted(set(jds))

#Loop over files and stokes 
prefix = args[0].split('.')[0]
for jd in jds:
    for s in stokes:
        file1 = prefix+'.'+jd+'.'+POL_WGTS[s].keys()[0]+'.'+file_ending
        file2 = prefix+'.'+jd+'.'+POL_WGTS[s].keys()[1]+'.'+file_ending
        outfile = prefix+'.'+jd+'.'+s+'.'+file_ending
        print file1, 'and',file2,'->', outfile
        if os.path.exists(outfile):
            print '    File exists, skipping.'
            continue
        if os.path.exists(file1) == False or os.path.exists(file2) == False:
            continue
        uvi1 = a.miriad.UV(file1)
        t2,d2,f2 = C.arp.get_dict_of_uv_data([file2],antstr='cross',polstr=POL_WGTS[s].keys()[1])
        
        uvi2 = a.miriad.UV(file2)
        uvo = a.miriad.UV(outfile, status='new')
        uvo.init_from_uv(uvi1, override={'pol':a.miriad.str2pol[s]})
        uvo.pipe(uvi1, mfunc=mfunc, raw=True, append2hist='COMBINE_POL:' + ' '.join(sys.argv) + '\n')
        try: uvo['var'] = uvi1['var']**2 + uvi2['var']**2 #if LST-binned already
        except: pass
        del(uvi1); del(uvo); del(uvi2)
        
