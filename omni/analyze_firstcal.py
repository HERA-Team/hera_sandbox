#! /usr/bin/env python
"""
Picking-out 'bad antennae' using a simple recursive 2sigma-clipping method, where
sigma is the standard deviation of the firstcal delay solutions over a night.
"""
import sys, optparse, numpy as np, aipy, glob
from matplotlib import pyplot as plt
#plt.use('Agg')

o = optparse.OptionParser()
o.add_option('--pcnt',default=93,help='Percentile cut to make on standard deviants per antenna delay. Default is worst 7% of antennae.')
o.add_option('--statfile',help='Full path and name of file to output mean and s.d. of delay solutions per antenna')
o.add_option('--badantsin',help='Full path and name of badantfile for previous analyze_firstcal run')
o.add_option('--badantsout',help='Full path and name of file to output list of bad antennae to')
#o.add_option('--makeplots',help='Full path to file to plot stats in.')
o.add_option('--verbose','-V',action='store_true',help='toggle verbosity')

opts,args = o.parse_args(sys.argv[1:])

#get initial parameters: assumes all npz files given have same badants
d = np.load(args[0])
DD = {}
for k in d.keys():
    if k.isdigit():
        DD[k] = []
jd = args[0].split('/')[-1].split('.')[1]
if opts.verbose: print jd
del(d)

for npzfile in args:
    print '    Reading %s'%npzfile
    d = np.load(npzfile)
    for k in d.keys():
        if k.isdigit():
            try: DD[k].append(d[k+'d'][0])
            except KeyError: 
                print 'Key error %s'%k
                continue

stds = []
#save antenna stats
if not opts.statfile is None: fstat = open(opts.statfile, 'w')
for k in DD.keys(): 
    stds.append(np.std(DD[k]))
    if not opts.statfile is None: print >> fstat, '%s %f %f'%(k, np.mean(DD[k]), np.std(DD[k]))
if not opts.statfile is None: fstat.close()

w = np.where(stds>=np.percentile(stds,float(opts.pcnt)))[0] #XXX best stat to cut on?
badants = list(np.array(DD.keys())[list(w)])
bastring=','.join(sorted(badants))

#doing it this way makes it more flexible: can run on a list of JDs, or just the one
if not opts.badantsin is None:
    if opts.verbose: print 'Using %s for input bad antennae'%opts.badantsin
    with open(opts.badantsin) as search:
        for line in search:
            if line.startswith(jd):
                   ba_1 = line.split(' ')[1]
                   if opts.verbose: print 'Found these bad antennae: %s'%ba_1
    bastring+=','+ba_1
    
if not opts.badantsout is None:
    print '    Saving badants to %s'%opts.badantsout
    fbad = open(opts.badantsout, 'w')    
    print >> fbad, '%s %s'%(jd,bastring)
    fbad.close()
else:
    if opts.verbose: print badants
