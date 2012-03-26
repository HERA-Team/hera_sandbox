#!/usr/bin/env python
#
#  fitmdl_mpi.py
#  
#
#  Created by Danny Jacobs on 3/24/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m,os
import sys, optparse,pickle,time,curses as C,logging,warnings
from mpi4py import MPI

"""
A 40% new implimentation of fitmdl using Message Passing Interface (MPI)
D. Jacobs
Averages over time.

Run with 
mpiexec -n <number of nodes> fitmdl_mpi.py [switches] <files> 
"""
o = optparse.OptionParser()
o.set_usage('fitmdl_mpi.py [options] *.uv')
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
o.add_option('-v',dest='verb',action='store_true',
    help="Print stuff.")
o.add_option('--vv',dest='vverb',action='store_true',
    help="Print even more")
opts, args = o.parse_args(sys.argv[1:])

#Logging bits
if opts.vverb: logging.disable(logging.DEBUG)
elif opts.verb: logging.disable(logging.INFO)
else: 
    logging.disable(logging.ERROR)
    warnings.simplefilter('ignore',Warning)
fname = sys.argv[0].split(os.sep)[-1]
log = logging.getLogger(fname)

# Parse command-line options
opts.ant += ',cross'
if opts.maxiter < 0: opts.maxiter = n.Inf

"basic fitmdl (no options)"
"features: stream scattering"
"send the array and catalog parameters to everyone"
"fit averaging overall deps (t,bl,freqs)"
#Begin: Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#End: Setup MPI


#setup curses
#myscreen = C.initscr()

first_fit=None 
def fit_func(prms):
    tstart = time.time()
    global first_fit
    if first_fit is None: print "first fit!"
    comm.bcast(1,root=0)
    comm.bcast(key_list,root=0)
    comm.bcast(prms,root=0)
    scores = comm.gather(0)
    m1s = comm.gather(0)
    score = n.sqrt(n.sum(scores)/n.sum(m1s))
    prms = a.fit.reconstruct_prms(prms, key_list)
    if not opts.quiet: a.fit.print_params(prms)
    if first_fit is None: 
        first_fit = score
        return score
    if not opts.quiet:
        print
        print 'Score:', n.real(score), 
        print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
        print '-' * 70
#        print "number of scores = %d, score = %f, m1 = %f, %f%%"%\
#            (len(scores)-1,n.sum(scores)/n.sum(m1s),n.sum(m1s),score/first_fit)
#        print "ttime = %5.2f"%(time.time()-tstart)
    return score/first_fit
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
    for (uvw,t,bl),d in uv.all():   
        if not bl in bls: bls.append(bl)
        if not t in ts: ts.append(t)
        nrec +=1
    ts = list(set(ts))
    bls = list(set(bls))
    uv.rewind()
    return ts,bls,freqs,nrec

#Begin: Setup environment

uv = a.miriad.UV(args[0])

#inshprms = None
#inprms = "cyg=jys"
#srcs = ['cyg']
#cats = ['misc']
#    print "opening ",'zen.2454565.41207.uvc'
#    uv = a.miriad.UV('zen.2454565.41207.uvc')
nchan = uv['nchan']

aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                     ('t', n.float), 
                     ('bl', '<i4', 2),
                     ('m1',n.complexfloating),
                     ('pol',str,2)]),
                ('spec', n.complexfloating,len(chans)),
                ('mask',n.bool,len(chans))])
#End: Setup Environment

if rank==0:  #The "master node"
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
    for prm in prms: prm_dict[prm] = prms[prm].keys()
    start_prms = aa.get_params(prm_dict)
    start_prms.update(cat.get_params(prm_dict))
    for obj in start_prms:
        for prm in start_prms[obj]:
            if prms[obj][prm][0] != None:
                start_prms[obj][prm] = prms[obj][prm][0]
            
    prm_list, key_list = a.fit.flatten_prms(start_prms)
    prm_list = n.array(prm_list)
    #End: Parse Parameters
    


    #divide along t only
    ts,bls,freqs,nrec = about(uv)
    tset = n.array_split(ts,size-1)
    #stream the data to the nodes
    print "streaming out data "
    print "tstart \t %15.4f"%(n.min(ts))
    print "tstop \t %15.4f"%n.max(ts)
    print "nbaselines %d"%len(bls)
    print "freqs (%4.3f : %4.3e : %4.3f)"%(n.min(freqs),n.diff(freqs)[0],n.max(freqs)) 
    print "nspectra %d"%nrec
    print "total data points %d"%(nrec*len(freqs))
    #this last is estimated assuming 2*8 bytes per complexfloat. Does not include preambles or flags.
    print "estimated total memory %5.2f"%(nrec*len(freqs)*len(args)*2*4.0/1024**3)
    print "estimated memory/process %5.2f"%(nrec*len(freqs)*len(args)*2*4.0/1024**3/(size-1))
    samp,vsamp,wgts = (0,0,0)
    for file in args:
        uv = a.miriad.UV(file)
        print "streaming: ", file
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        uv.select('decimate', opts.decimate, opts.decphs)
        for p,d,f in uv.all(raw=True):
            uvw,t,(i,j) = p
            if i==j: continue
            d = d.take(chans)
            f = f.take(chans)
            m2 = n.where(f,0,n.abs(d)**2).sum()
            p += (m2,)
            p += (a.miriad.pol2str[uv['pol']],)
            rec = n.array((p,d,f),dtype=uv_aipy,ndmin=1)
            #metrology
            samp += len(d)
            vsamp += n.logical_not(f).astype(n.int).sum()
            wgts += m2
            for i in range(1,size):
                if rec['preamble']['t'] in tset[i-1]:
                    comm.bcast(i,root=0)               #alert the grid about incoming
                    comm.send(rec, dest=i, tag=1)      #send the record
    comm.bcast(0,root=0)                           #alert the grid: finished buffering
    if not opts.quiet:
        print 'Cache summary:'
        print '   %d samples' % samp
        print '   %d valid' % vsamp
        print '   %f sum weights' %  wgts
    #set up the AntennaArray and Catalog
    uv.rewind()
#    cat_head = {'srcs':['cyg'],'catalogs':['misc']}
#    comm.bcast(cat_head,root=0)
    rv = a.optimize.fmin(fit_func,prm_list,
        full_output=1,disp=0,
        maxfun=opts.maxiter,maxiter=n.Inf,
        ftol=opts.ftol,xtol=opts.xtol)
    prms,score = rv[:2]
    prms = a.fit.reconstruct_prms(prms, key_list)
    print
    a.fit.print_params(prms)
    print 'Score:', score * first_fit, 
    print '(%2.2f%% of %f)' % (100 * score, first_fit)
    print '------------------------------------------------------------'
    comm.bcast(0,root=0)#clear the nodes
else: #All other nodes
    del(uv)
    data = n.array([],dtype=uv_aipy,ndmin=1)      #Create the data container
    #Buffer loop
    while(True):
        inc_rank = comm.bcast(root=0)             #Recieve instructions
        if inc_rank==0: break                     #Quit if we're done
        elif inc_rank==rank:                      #otherwise, recieve data, if its for us
            _in = comm.recv(source=0,tag=1)
            data = n.concatenate((data,_in))
    if len(data): print rank,"received ",len(data),"records having",len(data)*len(data[0]['spec']),"samples"
    else: print rank,"received no data"
#    cat_head = comm.bcast(None)
#    cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
    # sim difference loop
    print rank,"starting sim loop"
    while(True):
        if comm.bcast(root=0)==0: break
        tstart = time.time()
        key_list = comm.bcast(None,root=0)
        prm_list = comm.bcast(None,root=0)
        prms = a.fit.reconstruct_prms(prm_list, key_list)
        aa.set_params(prms)
        cat.set_params(prms)
        dra,ddec = cat.get('ionref')
        a1,a2,th = cat.get('srcshape')
        score = 0
        m1 = 0
        tnow = 0
        for p,d,f in data:
            if tnow != p['t']:
                tnow = p['t']
                aa.set_jultime(tnow)
                cat.compute(aa)
                aa.sim_cache(cat.get_crds('eq',ncrd=3),cat.get_jys(),
                    mfreqs=cat.get('mfreq'),
                    ionrefs=(dra,ddec),srcshapes=(a1,a2,th))
            i,j = p['bl']
#                print '\t',i,j
            if i==j: continue
            simd = aa.sim(i,j,pol=p['pol'])
            difsq = n.abs(d - simd)**2
            difsq = n.where(f,0,difsq)
            score += difsq.sum()
            m1 += p['m1']
        comm.gather(score)
        comm.gather(m1)
#        print "node %d, ptime = %5.2f"%(rank,time.time()-tstart)