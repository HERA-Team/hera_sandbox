#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys, string

o = optparse.OptionParser()
o.set_usage('omni_run.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--calpar',dest='calpar',type='string',
            help='Path and name of POL.p ("xx.p") calpar file from firstcal v1.')
o.add_option('--redinfo',dest='redinfo',type='string',default='',
            help='Path and name of .bin redundant info file from firstcal v1. If not given, it will be created.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .npz files. Include final / in path.')
o.add_option('--ba',dest='ba',default=None,
            help='Antennas to exclude, separated by commas.')
o.add_option('--fc2',dest='fc2', type='string',
            help='Path and name of POL.npz file outputted by firstcal v2 (comma delimited string if more than one polarization).')
o.add_option('--remove_degen',dest='rmdegen',default=True, 
            help='Toggle degeneracy removal functionality (recommended; default=True).')
o.add_option('--minV',action='store_true',
            help='Toggle V minimization capability. This only makes sense in the case of 4-pol cal, which will set crosspols (xy & yx) equal to each other')
o.add_option('--verbose',action='store_true',help='Toggle verbosity')
opts,args = o.parse_args(sys.argv[1:])

#Dictionary of calpar gains and files
pols = opts.pol.split(',')
files = {}
g0 = {} #firstcal gains

if opts.minV and len(list(set(''.join(pols))))==1:
    raise AssertionError('Stokes V minimization requires crosspols in the "-p" option.')

try: fc2 = opts.fc2.split(',')
except AttributeError: print('No fc2 files supplied. Continuing.')

for pp,p in enumerate(pols):
    #dictionary of calpars per pol
    if len(list(set(p))) > 1:
        if opts.verbose: print 'Not trying to get fcal info for pol %s'%p
        continue 
    g0[p[0]] = {} #indexing by one pol letter instead of two
    
    if opts.calpar != None: #if calpar is given
        if p in opts.calpar: #your supplied calpar matches a pol
            print 'Reading', opts.calpar
            cp = pickle.load(open(opts.calpar,'rb'))
            for i in xrange(cp[p].shape[1]): #loop over antennas
                g0[p[0]][i] = numpy.conj(cp[p][:,i] / numpy.abs(cp[p][:,i]))
        else: #looks for a calpar you haven't stated in the call
            new_cp = opts.calpar.split('.p')[0][:-2]+p+'.p' #XXX assumes calpar naming is *pol.p
            if os.path.exists(new_cp): #if it exists, use it
                print 'Reading', new_cp
                cp = pickle.load(open(new_cp,'rb'))
                for i in xrange(cp[p].shape[1]): #loop over antennas
                    g0[p[0]][i] = numpy.conj(cp[p][:,i] / numpy.abs(cp[p][:,i]))
            elif len(list(set(p))) > 1: #if the crosspol first_cal is missing, don't worry
                #print '%s not found, but that is OK'%new_cp
                continue
            else: #if the linpol first_cal is missing, do worry
                raise IOError('Missing first_cal file %s'%new_cp)
    if opts.fc2: #if fc2 file is given
        fc2file = next((s for s in fc2 if p in s), None)
        if not fc2file == None:
            print 'Reading %s, pol=%s'%(fc2file,p)
            _,_g0,_,_ = capo.omni.from_npz(fc2file)
            for i in _g0[p[0]].keys():
                g0[p[0]][i] = _g0[p[0]][i][:,:] #/ numpy.abs(_g0[p[0]][i][:,:])
        elif len(list(set(p))) > 1:
            continue #don't use crosspols to firstcal
        else:
            raise IOError("Please provide a valid first cal file for polarization %s"%p) 

for filename in args:
    files[filename] = {}
    for p in pols:
        fn = filename.split('.')
        fn[3] = p
        files[filename][p] = '.'.join(fn)

#Create info
if opts.redinfo != '': #reading redinfo file
    print 'Reading',opts.redinfo
    info = omnical.info.RedundantInfoLegacy()
    if opts.verbose: print '   Getting reds from redundantinfo'
    info.fromfile(opts.redinfo)
else: #generate reds from calfile
    aa = aipy.cal.get_aa(opts.cal,numpy.array([.15]))
    if opts.verbose: print 'Getting reds from calfile'
    if opts.ba: #XXX assumes exclusion of the same antennas for every pol
        ex_ants = []
        for a in opts.ba.split(','):
            ex_ants.append(int(a))
        print 'Excluding antennas:',sorted(ex_ants)
    else: ex_ants = []
    info = capo.omni.aa_to_info(aa, pols=list(set(''.join(pols))), ex_ants=ex_ants, crosspols=pols, minV=opts.minV)
reds = info.get_reds()

### Omnical-ing! Loop Through Compressed Files ###
for f,filename in enumerate(args):
    file_group = files[filename] #dictionary with pol indexed files
    print 'Reading:'
    for key in file_group.keys(): print '   '+file_group[key]
    if len(pols)>1: #zen.jd.npz
        npzb = 3
    else: #zen.jd.pol.npz
        npzb = 4 
    npzname = opts.omnipath+'.'.join(filename.split('/')[-1].split('.')[0:npzb])+'.npz'
    if os.path.exists(npzname):
        print '   %s exists. Skipping...' % npzname
        continue

    timeinfo,d,f = capo.arp.get_dict_of_uv_data([file_group[key] for key in file_group.keys()], antstr='cross', polstr=opts.pol)
    t_jd = timeinfo['times']
    t_lst = timeinfo['lsts']
    freqs = numpy.arange(.1,.2,.1/len(d[d.keys()[0]][pols[0]][0]))
    SH = d.values()[0].values()[0].shape #shape of file data (ex: (19,203))
    data,wgts,xtalk = {}, {}, {}
    m2,g2,v2 = {}, {}, {}
    for p in g0.keys():
        for i in g0[p]: 
            if g0[p][i].shape != (len(t_jd),len(freqs)):
                g0[p][i] = numpy.resize(g0[p][i],SH) #resize gains like data
            else: continue
    data = d #indexed by bl and then pol (backwards from everything else)
    for p in pols:
        wgts[p] = {} #weights dictionary by pol
        for bl in f: 
            i,j = bl
            wgts[p][(j,i)] = wgts[p][(i,j)] = numpy.logical_not(f[bl][p]).astype(numpy.int)
    if len(g0[g0.keys()[0]]) == 0:
        if opts.verbose: print "   Making initial gains of all 1's..."
        for key in g0:
            g0[key] = {}
            for ant in info.subsetant:
                g0[key][ant] = numpy.ones(SH,dtype='complex')
    print '   Logcal-ing' 
    m1,g1,v1 = capo.omni.redcal(data, info, gains=g0, removedegen=False)
    print '   Lincal-ing'
    m2,g2,v2 = capo.omni.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=True)
    xtalk = capo.omni.compute_xtalk(m2['res'], wgts) #xtalk is time-average of residual
    m2['history'] = 'OMNI_RUN: '+' '.join(sys.argv) + '\n'
    m2['jds'] = t_jd
    m2['lsts'] = t_lst
    m2['freqs'] = freqs
    
    print '   Saving %s'%npzname
    capo.omni.to_npz(npzname, m2, g2, v2, xtalk)
    
