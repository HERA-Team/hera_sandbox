#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys
import numpy as np

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply_fc.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,pol=True)
o.add_option('--suffix',dest='suffix',type='string',default='',help='suffix for fcal directory, e.g. _allEW')
o.add_option('--round', dest='round', type='int',
            help='Round of firstcal flagging used. Assumes filestructure /path/to/data/fcal_{pol}_{round}_{suffix}/*fc.npz')
opts,args = o.parse_args(sys.argv[1:])

pols = opts.pol.split(',')
assert(len(pols)==1) #this script is for PAPER-128 fcs, which are one per file. We're not swapping-out pols

### Read Data and Solutions ###
for i,filename in enumerate(args):
        if not len(opts.suffix)>0:
            fcfile = os.path.dirname(filename)+'/fcal_%s_%i/%s.fc.npz'%(opts.pol,opts.round,os.path.basename(filename))
        else:
            fcfile = os.path.dirname(filename)+'/fcal_%s_%i%s/%s.fc.npz'%(opts.pol,opts.round,opts.suffix,os.path.basename(filename))
        
        if not os.path.exists(fcfile):
            print 'No firstcal file %s ... Skipping...'%fcfile
            continue
        
        fcdata = np.load(fcfile)
        gains = {}
        for k in fcdata.keys():
            if k.endswith(opts.pol[0]):
                gains[int(k[:-1])] = fcdata[k][0,:]
        
        newfile = filename+'f' #saves in same place as uvcRRE file
        if os.path.exists(newfile):
            print '    %s exists.  Skipping...' % newfile
            continue
        
        times = []
        def mfunc(uv,p,d,f): #loops over time and baseline
            global times #global list
            _,t,(a1,a2) = p #parse preamble
            if len(times) == 0 or times[-1] != t: times.append(t) #fill times list
            ti = len(times) - 1 #time index
            try: d*=numpy.conj(gains[a1])
            except(KeyError): pass
            try: d*=numpy.conj(gains[a2].conj())
            except(KeyError): pass
            return p, numpy.where(f,0,d), f
        
        print '    Applying fc solns to %s'%filename 
        uvi = aipy.miriad.UV(filename)
        uvo = aipy.miriad.UV(newfile,status='new')
        uvo.init_from_uv(uvi)
        print '    Saving', newfile
        uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='omni_apply_fc: ' + ' '.join(sys.argv) + '\n')
        
