#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,pol=True)
o.add_option('--xtalk',dest='xtalk',default=False,action='store_true',
            help='Toggle: apply xtalk solutions to data. Default=False')
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
o.add_option('--firstcal', action='store_true', 
            help='Applying firstcal solutions.')
o.add_option('--ubls', default=[],
            help='List of unique baselines to include in *O file, separated by semi-colons (ex: "0,2;0,1"). The default will include all baselines.')
o.add_option('--ba',dest='ba',default=None,
            help='Antennas to exclude, separated by commas.')
opts,args = o.parse_args(sys.argv[1:])

#File Dictionary
pols = opts.pol.split(',')
files = {}
for filename in args:
    files[filename] = {}
    for p in pols:
        fn = filename.split('.')
        fn[3] = p
        files[filename][p] = '.'.join(fn)

#List of UBLS and ANTS
aa = aipy.cal.get_aa('psa6622_v003',numpy.array([.15])) #XXX PSA128
_,_,bl2sep = capo.zsa.grid2ij(aa.ant_layout)
if opts.ubls != []:
    keep = opts.ubls.split(';')
else: keep = []

### Read Data and Solutions ###
for f,filename in enumerate(args):
    if len(pols)>1:
        npzb=3
    else:
        npzb=4
    omnifile = opts.omnipath % '.'.join(filename.split('/')[-1].split('.')[0:npzb])
    print '   Omnical npz:', omnifile
    _,gains,_,xtalk = capo.omni.from_npz(omnifile) #loads npz outputs from omni_run
    for p in pols:
        print 'Reading', files[filename][p]
        if opts.firstcal:
            newfile = files[filename][p].split('/')[-1]+'F' #saves in cwd. Firstcal ext.
        else:
            newfile = '/'.join(omnifile.split('/')[:-1])+'/'+files[filename][p].split('/')[-1]+'O' #saves wherever the omnifile is
        if os.path.exists(newfile):
            print '    %s exists.  Skipping...' % newfile
            continue
        times = []
    
        def mfunc(uv,p,d,f): #loops over time and baseline
            global times #global list
            _,t,(a1,a2) = p
            p1,p2 = pol = aipy.miriad.pol2str[uv['pol']]
            if a1==a2: return p,None,None #skip autos
            try: trysep = bl2sep[aipy.miriad.ij2bl(a1,a2)]
            except: return p,None,None #outriggers
            if trysep not in keep and keep != []: return p,None,None #skip some baselines if specified
            if opts.ba != None:
                if a1 in map(int,opts.ba.split(',')) or a2 in map(int,opts.ba.split(',')): return p,None,None #skip some antennas if specified
            if len(times) == 0 or times[-1] != t: times.append(t) #fill times list
            if opts.xtalk: #subtract xtalk
                try:
                    xt = numpy.resize(xtalk[pol][(a1,a2)],d.shape)
                    d -= xt
                except(KeyError):
                    try:
                        xt = numpy.resize(xtalk[pol][(a2,a1)].conj(),d.shape) 
                        d -= xt
                    except(KeyError): pass
            if opts.firstcal:
                ti = 0
            else:
                ti = len(times) - 1 #time index
            try: d /= gains[p1][a1][ti] #apply gains
            except(KeyError): pass
            try: d /= gains[p2][a2][ti].conj() 
            except(KeyError): pass
            return p, numpy.where(f,0,d), f
    
        if opts.xtalk: print '    Calibrating and subtracting xtalk'
        else: print '    Calibrating'
        uvi = aipy.miriad.UV(files[filename][p])
        uvo = aipy.miriad.UV(newfile,status='new')
        uvo.init_from_uv(uvi)
        print '    Saving', newfile
        uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='OMNICAL: ' + ' '.join(sys.argv) + '\n')
        
