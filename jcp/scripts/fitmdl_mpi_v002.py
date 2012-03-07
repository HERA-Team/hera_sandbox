#!/usr/bin/env python
#
#  fitmdl_mpi_v002.py
#  
#
#  Created by Danny Jacobs on 8/15/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pickle,time,re,glob,os,datetime
import mpi4py.rc,sys,logging,warnings
from copy import copy,deepcopy
mpi4py.rc.threaded = False
from mpi4py import MPI
tstart = time.time()
"""
A 40% new implimentation of fitmdl using Message Passing Interface (MPI)
D. Jacobs
Averages over time.

Run with 
mpiexec -n <number of nodes> fitmdl_mpi.py [switches] <files> 
"""
o = optparse.OptionParser()
o.set_usage('beamform_mpi_test.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True,
    cal=True, src=True, dec=True, prms=True)
o.add_option('-S', '--shared_prms', dest='shprms',
    help='Parameter listing w/ same syntax as "-P/--prms" except that all objects listed for a parameter will share an instance of that parameter.')
o.add_option('-q', '--quiet', dest='quiet', action='store_true',
    help='Be less verbose.')
o.add_option('--maxiter', dest='maxiter', type='float', default=-1,
    help='Maximum # of iterations to run.  Default is infinite.')
o.add_option('--xtol', dest='xtol', type='float', default=1e-10,
    help='Fractional change sought in it parameters before convergence.  Default 1e-10.')
o.add_option('--ftol', dest='ftol', type='float', default=1e-10,
    help='Fractional tolerance sought in score before convergence.  Default 1e-10.')
o.add_option('-b', '--beam', dest='beam', action='store_true',
    help='Normalize by the primary beam response in the direction of the specified source.')
o.add_option('-f', '--srcflux', dest='srcflux', action='store_true',
    help='Normalize by the spectrum of the specified source.')
o.add_option('--minuv', dest='minuv', type='float', default=0,
    help='Minimum uv length (in wavelengths) for a baseline to be included.')    
#o.add_option('--gsm',
#    help="Fit all modelled sources. Specify source parameters to fit. comma list of any (jys,dra,ddec,a1,a2,th)")
#o.add_option('--tracefile',
#    help="Output a text file with snapshot parameters, only works with --gsm option. Sources are demarkaded by \\t and variables by ,")
#o.add_option('--filter',action='store_true',
#    help="Filter the data before doing the beamform.")
#o.add_option('-d', '--dw', dest='dw', type=int, 
#    help='The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
o.add_option('--name',default='beamform_mpi',
    help="specify an experiment name, will output trace-files for all parameters.")
o.add_option('--time',action='store_true',
    help="Log code execution times. Requires an experiment --name")
o.add_option('--dt',type='str',default='1n',
    help='Amount of time to step at a time. Give number with units (n or s eg integration number or seconds). default 1n')
o.add_option('-v',dest='verb',action='store_true',
    help="Print stuff.")
o.add_option('--vv',dest='vverb',action='store_true',
    help="Print even more")

opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
opts.ant += ',cross'
if opts.maxiter < 0: opts.maxiter = n.Inf
if opts.dt.endswith('n'):
    dt = int(opts.dt[:-1])
else:
    dt = float(opts.dt[:-1])

if opts.vverb: logging.disable(logging.DEBUG)
elif opts.verb: logging.disable(logging.INFO)
else: 
    logging.disable(logging.ERROR)
    warnings.simplefilter('ignore',Warning)
fname = sys.argv[0].split(os.sep)[-1]
log = logging.getLogger(fname)
    

#Begin: Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#End: Setup MPI

first_fit=None 
def fit_func(prms):
    global first_fit
    if first_fit is None: print "first fit!"
    comm.bcast(1,root=0)
    comm.bcast(key_list,root=0)
    comm.bcast(prms,root=0)
    scores = comm.gather(0)
    m1s = comm.gather(0)
#    bmflux = n.sum(comm.gather(0))
#    print "bmflux = ",bmflux
    score = n.sqrt(n.sum(scores)/n.sum(m1s))
    prms = a.fit.reconstruct_prms(prms, key_list)
    if not opts.quiet: a.fit.print_params(prms)
#    if not opts.quiet: print "number of scores = %d, score = %f, m1 = %f, %%=%f"%(len(scores)-1,n.sum(scores)/n.sum(m1s),n.sum(m1s),score/first_fit)
    if first_fit is None: first_fit = score
    cat.set_params(prms)
    pprms = cat.get_params()
    fitlog([pprms,score],t,uv['lst'],final=False)
    if not opts.quiet:
        print
        print 'Score:', score, 
        print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
        print '-' * 70
    return score
def about(uv):
    """
    Read a miriad.UV and return a list of times, baselines, frequencies and total sample count
    Plays nice and rewinds uv when finished
    """
    bls = []
    ts = []
    nchan = uv['nchan']
    sfreq = uv['sfreq']
    dfreq = uv['sdf']
    freqs = n.arange(uv['sfreq'],uv['sfreq']+uv['sdf']*uv['nchan'],uv['sdf'])
    nrec=0.0
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    for (uvw,t,bl),d in uv.all():   
        if not bl in bls: bls.append(bl)
        if not t in ts: ts.append(t)
        nrec +=1
    ts = list(set(ts))
    bls = list(set(bls))
    uv.rewind()
    return ts,bls,freqs,nrec
def beamform(srcs):
    comm.bcast(3,root=0)#trigger flux computation]
    comm.bcast(srcs,root=0)
    fluxes =  comm.gather(None,root=0)
    fluxes = n.ma.array(fluxes[1:])
    cnt = comm.gather(None,root=0)
    cnt = n.ma.array(cnt[1:])
    return n.ma.sum(fluxes,axis=0)/n.ma.sum(cnt.clip(1,n.Inf),axis=0)
def exp_number(match):
    "match=ls explicit search criteria, eg name_*.txt, where * is replaced with next highest number startinga at 1" 
    files = glob.glob(match)
    if len(files)<1:
        return 1
    numlist = []
    for f in files:
        fnude = os.path.splitext(f)[0]
        try: 
            num = re.findall('[0-9]+$',fnude)[0]
            numlist.append(int(num))
        except(IndexError): pass
    numlist.sort()
    if len(numlist)==0: return 1
    return numlist[-1]+1
def fitlog(prms_score,jd,lst,final=False,extra=[]):
    prms,score = prms_score
    pflat,kys = a.fit.flatten_prms(prms)
    if not os.path.exists(root+'/fit_kys.pkl'):
        pickle.dump(kys,open(root+'/fit_keys.pkl','w'))
#    open(root+'/fit_keys.txt','a').write(pickle.dumps(pflat)+'\n')
    open(root+'/fit_times.txt','a').write(str(jd)+'\t'+str(lst)+'\n')
    open(root+'/fit_scores.txt','a').write(str(int(final))+'\t'+str(score)+'\t'+str(first_fit)+'\n')
    open(root+'/fit_prms.txt','a').write('\t'.join(map(str,pflat))+'\n')
def bmlog(prms,jd,lst,extra=[]):
    pflat,kys = a.fit.flatten_prms(prms)
    if not os.path.exists(root+'/bm_kys.pkl'):
        pickle.dump(kys,open(root+'/bm_keys.pkl','w'))
    open(root+'/bm_keys.txt','a').write(pickle.dumps(pflat)+'\n')
    open(root+'/bm_times.txt','a').write(str(jd)+'\t'+str(lst)+'\n')
    open(root+'/bm_prms.txt','a').write('\t'.join(map(str,pflat))+'\n')
def timelog(file,t,name):
    if not os.path.exists(file):
        open(file,'w').write('#'+'\t'.join(['name','t'])+'\n')
    open(file,'a').write('\t'.join(map(str,[name,t]))+'\n')
def srclog(srcs):
    allsrcs = cat.keys()
    allsrcs.sort()
    if not os.path.exists(root+'/srcs.txt'):
        open(root+'/srcs.txt','w').write('\t'.join(allsrcs)+'\n') 
    T = []
    for src in allsrcs:
        if src in srcs: T.append('1')
        else: T.append('0')
    open(root+'/srcs.txt','a').write('\t'.join(T)+'\n')
#Begin: Setup environment
uv = a.miriad.UV(args[0])
nchan = uv['nchan']

aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
mfq = cat.get('mfreq')

uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                     ('t', n.float), 
                     ('bl', '<i4', 2),
                     ('m1',n.float),
                     ('pol',str,2)]),
                ('spec', n.complex64,len(chans)),
                ('mask',n.bool,len(chans))])
ts,bls,freqs,nrec = about(uv)

#tset = n.array_split(ts,size-1)
#End: Setup Environment
del(uv)   
#log my executions!


if rank==0:  #The "master node"
    print "streaming out data "
    print "divining stream parameters from ",args[0]
    if type(dt) is int:
        dt = float(dt)*n.diff(ts)[0]
    nrec *= len(args)
    print "averaging over %15.4f seconds"%dt



    print "tstart \t %15.4f"%(n.min(ts))
#    print "tstop \t %15.4f"%n.max(ts)
    
    print "nbaselines %d"%len(bls)
    print "freqs (%4.3f : %4.3e : %4.3f)"%(n.min(freqs),n.diff(freqs)[0],n.max(freqs)) 
    print "nspectra %d"%nrec
    print "total data points %d"%(nrec*len(freqs))
    #this last is estimated assuming 2*8 bytes per complexfloat. Does not include preambles or flags.
    print "estimated total memory %5.2f"%(nrec*len(freqs)*len(args)*2*4.0/1024**3)
    print "estimated memory/process %5.2f"%(nrec*len(freqs)*len(args)*2*4.0/1024**3/(size-1))
    #if the time was given as step number, find the number of times.
    
    
    #initialize trace files
    if not opts.name is None:
        #find the fit trace file name
        #get the experiment number for this name
        enum = exp_number(opts.name+'*')
        
        fitname = opts.name+'_fit_'+str(enum)+'.txt'
        #start a beam trace
        bmname = opts.name+'_bm_'+str(enum)+'.txt'
        print "beamform_mpi: logging fit trace in %s and beamform results in %s"%(fitname,bmname)
        if opts.time:
            timefile = opts.name +'_time_'+str(enum)+'.txt'
        root = opts.name+'_'+str(enum)
        if not os.path.exists(root):
            print "creating directory %s"%root
            os.mkdir(root)
    comment = "mpi job: %s,%i with %d nodes on %s"%(opts.name,enum,size,os.uname()[1]) 
    open(os.getcwd()+os.sep+sys.argv[0].split(os.sep)[-1]+'.log','a').write(time.strftime("%D %H:%M:%S") +\
    '\t'+ ' '.join(sys.argv)+\
    '\t'+comment+'\n')

    #Begin: Parse Parameters

    prms, prm_dict, shkeys = {}, {}, []
    # Handle shared parameters
    if opts.shprms:
        shprms = map(a.scripting.parse_prms, opts.shprms.split(','))
        for s in shprms:
            keys = s.keys(); keys.sort()
            k = keys[0]
            # Only enter one instance of shared parameters (indexed under first key)
            if prms.has_key(k): prms[k].update(s[k])
            else: prms.update({k:s[k]})
            # Add an entry to shkeys for propagating variables to other objects
            shkeys.append((keys, s[k].keys()))
    # Handle normal parameters
    if not opts.prms is None:
        pd = a.scripting.parse_prms(opts.prms)
        for k in pd:
            if prms.has_key(k): prms[k].update(pd[k])
            else: prms[k] = pd[k]
    for prm in prms: prm_dict[prm] = prms[prm].keys();print prm,prms[prm].keys()
#
#    for obj in start_prms:
#        for prm in start_prms[obj]:
#            if prms[obj][prm][0] != None:
#                start_prms[obj][prm] = prms[obj][prm][0]
#    gsm_dict = dict(zip([cat[s].src_name for s in cat],[opts.gsm.split(',')]*len(cat)))

    start_prms = aa.get_params(prm_dict)
    start_prms.update(cat.get_params(prm_dict))
    cur_prms = start_prms.copy()
    #End: Parse Parameters
    print start_prms

    #stream the data to the nodes
    print "streaming out data "
    print "fitting will occur as we work through the stream"
    print dt,nrec
    tnow=1e9
    first=True
    data = n.array([],dtype=uv_aipy,ndmin=1)
    scatter_time = True
    startup_time = True
    read_time = True
    fit_time = True
    loop_time = True
    tnext = n.min(ts) + dt
    files = n.sort(args)
    rnum = 0
    newt = True
    for file in files:
        uv = a.miriad.UV(file)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        uv.select('decimate', opts.decimate, opts.decphs)
        tread = time.time()
        for p,d in uv.all():
            #data prep and loop init stuff
            rnum +=1
            d = d.take(chans)
            m2 = n.float(n.where(d.mask,0,n.ma.abs(d)**2).sum())
            p += (m2,)
            p += (a.miriad.pol2str[uv['pol']],)
            t = p[1]
            if len(chans)!=1: 
                mask = d.mask
                d = d.data
            else:
                mask = d.mask[0]
                d = d.data[0]
                
            #select buffer/op operation. Optinos are 
            #       "very first record", 
            #        "within current buffered period"                            
            #        "last record file or 
            if t<tnext and not rnum==nrec:
                tnow = t
                lstnow = uv['lst']
                if newt:
#                    print "starting stream at time %15.4f"%t
                    print '   begin   %f, lst   %f'%(t,lstnow)
                    tstart = t
                    newt = False
                if opts.time and read_time: timelog(timefile,'rec_read',time.time()-tread); read_time=False;
                if opts.time and loop_time: tloop = time.time()
                rec = n.array((p,d,mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,rec))
            if (tnext == t or rnum==nrec):
                #if this is the last record, add it to the pile.
                if rnum==nrec:
                    rec = n.array((p,d,mask),dtype=uv_aipy,ndmin=1)
                    data = n.concatenate((data,rec))

                #scatter the data, taking care not to send any to rank=0
                #Obfuscated by a bunch of timing metrology.
                if opts.time and startup_time: 
                    timelog(timefile,'startup',time.time()-tstart)  #The amount of time to start the script and get to sending out the first data chunk
                    startup_time = False
                comm.bcast(2,root=0) #Send the "load data" command
                if opts.time and scatter_time:
                    tscatter = time.time()
                comm.scatter([0]*1+n.array_split(data,size-1),root=0)
                if opts.time and scatter_time:
                    timelog(timefile,'scatter',time.time()-tscatter) #the amount of time to scatter the data
                    scatter_time = False

                print "stream cache/interval summary"
                print '   end   %f, lst   %f'%(t,lstnow)
                print '   length %f seconds'%(t-tstart)
                print '   %d/%d'%(rnum,nrec)
                print '   %d samples'%n.sum(n.ones_like(data['spec']))
                print '   %d valid'%n.sum(n.logical_not(data['mask'].astype(n.int)))
                print '   %f sum weights'%n.sum(data['preamble']['m1'])
                #ok, we have our data. Now lets do something with it.

                #set up the instrument and sky models
                aa.set_jultime(p[1])
                cat.compute(aa)
                #reset the nodes with the complete catalog
                comm.bcast(4,root=0)
                comm.bcast(cat.keys(),root=0)                
                cur_prms.update(start_prms)
                srcs = [src for src in cat if cat[src].alt>0]
                prm_list, key_list = a.fit.flatten_prms(cur_prms)
                prm_list = n.array(prm_list)

                #do the fit for these times
                first_fit=None
                if opts.time and fit_time: tfit = time.time()
                rv = a.optimize.fmin(fit_func,prm_list,
                    full_output=1,disp=0,
                    maxfun=opts.maxiter,maxiter=n.Inf,
                    ftol=opts.ftol,xtol=opts.xtol)              
                if opts.time and fit_time: 
                    timelog(timefile,'complete fit',time.time()-tfit)
                    fit_time=False
                #mash up the results for printing and logging
                prms,score = rv[:2]
                prms = a.fit.reconstruct_prms(prms, key_list)

                fitlog([prms,score],tnow,lstnow,final=True)   #LOG EVENT
                a.fit.print_params(prms)
                print
                print 'Score:', score, 
                print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
                print '-' * 70
                tnext += dt
                lstnow=uv['lst']
                newt = True
                if opts.time and loop_time:
                    timelog(timefile,'loop',time.time()-tloop)
                    loop_time = False
                    
                #if we're to go on to a new time-step. Clean out our data buffer, put in the single t=tnext record and get ready.
                if tnext <= n.max(ts):
                    data = n.array([],dtype=uv_aipy,ndmin=1) #clear the time buffer
                    rec = n.array((p,d,mask),dtype=uv_aipy,ndmin=1)
                    data = n.concatenate((data,rec))
    comm.bcast(0,root=0) #reset the other nodes
    print "This was %s # %s in files %s"%(root,enum,' '.join(args))
else: #All other nodes
    print "node",rank,"starting up"
    data = None
    good_srcs = cat.keys()#a priori, they are all good!
    while(True):
        #get the latest command
        cmd = comm.bcast(root=0)
        if cmd==4:#update ok source list, and catalog
            good_srcs = comm.bcast(root=0)
            cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs) #restore catalog to pristine
            for src in cat.keys():#clean out bad sources
                if src not in good_srcs: 
                    print rank,"removing",cat.pop(src)
        if cmd==3: #Use current model to compute beamformed flux
            srcs = comm.bcast(root=0) #list of {name:RFB} dicts
            #beamform part
            fluxes = []
            cnts   = []
            Flags = []
            _aa = a.fit.AntennaArray((aa.lat,aa.long),aa.ants)
            ts_this_chunk = []
#            cat = ocat
            for snum,name in enumerate(srcs):
                tsim = 0
                for p,d,f in data:
                    t = p['t']
                    i,j = p['bl']
                    if t != tsim:
                        F = f
                        src = cat[name] 
                        aa.set_jultime(t)
                        cat.compute(aa)
                        #catch the sky model sim
                        flux = n.zeros_like(data['spec'][0]).astype(n.float)
                        cnt = n.zeros_like(flux).astype(n.int)
                        tmfq = cat.get('mfreq')
                        eqs = cat.get_crds('eq', ncrd=3)
                        flx = cat.get_jys()
                        dra,ddec = cat.get('ionref')
                        a1,a2,th = cat.get('srcshape')
                        aa.sim_cache(eqs, flx, mfreqs=tmfq, 
                                ionrefs=(dra,ddec), srcshapes=(a1,a2,th))
                        #cache the addition simulation
                        _aa.set_jultime(t)
                        _cat = a.fit.SrcCatalog([cat[name]])
                        _eqs = _cat.get_crds('eq', ncrd=3)
                        _flx = _cat.get_jys()
                        _dra,_ddec = _cat.get('ionref')
                        _a1,_a2,_th = _cat.get('srcshape')
                        _mfq = _cat.get('mfreq')
                        _aa.sim_cache(_eqs, _flx, mfreqs=_mfq, 
                                ionrefs=(_dra,_ddec), srcshapes=(_a1,_a2,_th))                    
                        tsim = t
                        if not t in ts_this_chunk: ts_this_chunk.append(t)
                    try: 
                        F = n.logical_or(F,f)
                        d = aa.phs2src(d, src, i, j)
                        if not opts.dw is None:#do delay transform
                            y1, y2 = opts.dw, -opts.dw
                            if y2 == 0: y2 = d.shape[1]
                            #subtract ->  delay filter -> beamform
                            d = n.where(f, 0, d)
                            flags = n.logical_not(f).astype(n.float)
                            ker = n.fft.ifft(flags)
                            sim_add = _aa.sim(i,j,pol=p['pol'])
                            sim_add = _aa.phs2src(sim_add, src, i, j)
                            sim_sub = aa.sim(i,j,pol=p['pol'])
                            sim_sub = aa.phs2src(sim_sub, src, i, j)
                            d = d-sim_sub #-subtract- divide by model
                            gain = n.sqrt(n.average(flags**2))
                            d = n.fft.ifft(d) 
                            if not n.all(d == 0):
                                d, info = a.deconv.clean(d, ker, tol=opts.clean)
                                d += info['res'] / gain
                            d[y1:y2] = 0
                            d = n.fft.fft(d)+sim_add
                        u,v,w = aa.gen_uvw(i, j, src)
                        tooshort = n.where(n.sqrt(u**2+v**2) < opts.minuv, 1, 0).squeeze()
                        if n.all(tooshort):
                            #print i,j, 'too short:',
                            #print n.average(n.sqrt(u**2+v**2)), '<', opts.minuv
                            continue                        
                        f = n.logical_or(f, tooshort)
                        gain = _aa.passband(i,j)
                        if opts.beam: gain *= _aa.bm_response(i,j,pol=opts.pol).squeeze()
                        d /= gain                       
                    except(a.phs.PointingError): 
                        d *= 0
                        print rank,"%s not up at %f"%(name,p['t'])#len(aa.bm_response(i,j,pol=opts.pol)[:,0]));d *=0
                    flux += d
                    cnt += n.logical_not(f).astype(n.int)
                fluxes.append(flux)
                cnts.append(cnt.clip(1,n.Inf))
                Flags.append(F)
            comm.gather(n.ma.array(fluxes,ndmin=2,mask=Flags),root=0)
            comm.gather(n.ma.array(cnts,ndmin=2,mask=Flags)*len(ts_this_chunk),root=0)                    
        elif cmd==2: 
            data=comm.scatter(root=0)
            continue 
            """Compute a score on currently buffered data. Waits for broadcast of 
            keys and parameters as output by a.fit.flatten_params().
            Returns a score and the "weights" via comm.gather in that order.
            """
        elif cmd==1:
            key_list = comm.bcast(None,root=0)
            prm_list = comm.bcast(None,root=0)
            #------
            
            prms = a.fit.reconstruct_prms(prm_list, key_list)
    #        a.fit.print_params(prms)
            try: aa.set_params(prms)
            except(IndexError):pass
            try: cat.set_params(prms)
            except(IndexError):pass
            dra,ddec = cat.get('ionref')
            a1,a2,th = cat.get('srcshape')
            score = 0
            m1 = 0
            tnow = 0
            for p,d,f in data:
                i,j = p['bl']
                if tnow != p['t']:
                    tnow = p['t']
                    aa.set_jultime(tnow)
                    cat.compute(aa)
                    eqs = cat.get_crds('eq', ncrd=3)
                    flx = cat.get_jys()
                    dra,ddec = cat.get('ionref')
                    aa.sim_cache(eqs, flx, mfreqs=mfq, 
                        ionrefs=(dra,ddec), srcshapes=(a1,a2,th))
                if i==j: continue
                simd = aa.sim(i,j,pol=p['pol'])
                difsq = n.abs(d - simd)**2
                difsq = n.where(f,0,difsq)
                score += difsq.sum()
                m1 += p['m1']
            comm.gather(score)
            comm.gather(m1)
        elif cmd==0: sys.exit() #exit loop