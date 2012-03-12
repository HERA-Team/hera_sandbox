#!/usr/bin/env python
#
#  fitmdl_mpi_snap.py
#  
#
#  Created by Danny Jacobs on 3/25/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pickle,time
from mpi4py import MPI

"""
A 40% new implimentation of fitmdl using Message Passing Interface (MPI).
D. Jacobs
Fits as a function of time (eg snap mode)


Run with 
mpiexec -n <number of nodes> fitmdl_mpi.py [switches] <files> 
"""
o = optparse.OptionParser()
o.set_usage('fitmdl_mpi_snap.py [options] *.uv')
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
o.add_option('--gsm',
    help="Fit all modelled sources. Specify source parameters to fit. comma list of any (jys,dra,ddec,a1,a2,th)")
o.add_option('--tracefile',
    help="Output a text file with snapshot parameters, only works with --gsm option. Sources are demarkaded by \\t and variables by ,")
    
opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
opts.ant += ',cross'
if opts.maxiter < 0: opts.maxiter = n.Inf


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
    score = n.sqrt(n.sum(scores)/n.sum(m1s))
    prms = a.fit.reconstruct_prms(prms, key_list)
    if not opts.quiet: a.fit.print_params(prms)
#    if not opts.quiet: print "number of scores = %d, score = %f, m1 = %f, %%=%f"%(len(scores)-1,n.sum(scores)/n.sum(m1s),n.sum(m1s),score/first_fit)
    if first_fit is None: first_fit = score
    if not opts.quiet:
        print
        print 'Score:', score, 
        print '(%2.2f%% of %f)' % (100 * score / first_fit, first_fit)
        print '-' * 70
    return score / first_fit
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
def print_trace(rv,prm_dict):
    """
    Prints a line source parameters from a optimize return value.
    Pairs nicely with opt.gsm and opts.tracefile
    """
    global first_fit
    prms,score = rv[:2]
    prms = a.fit.reconstruct_prms(prms, key_list)
#    print prms,prm_dict
    out = ''
    for name,parms in prm_dict.iteritems():
        for parm in parms: out += str(prms[name][parm])+'\t'
        out += "\t"
    out += str(score*first_fit)+"\t"+str(100*score)+"\t" +str(first_fit) +"\n"
    return out

    
    
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

tset = n.array_split(ts,size-1)
#End: Setup Environment
del(uv)
if rank==0:  #The "master node"
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

    #Begin: Parse Parameters
    if opts.gsm is None:
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
    else:
        prm_dict = dict(zip([cat[s].src_name for s in cat],[opts.gsm.split(',')]*len(cat)))
    start_prms = aa.get_params(prm_dict)
    start_prms.update(cat.get_params(prm_dict))
#    print start_prms
#    for obj in start_prms:
#        for prm in start_prms[obj]:
#            if prms[obj][prm][0] != None:
#                start_prms[obj][prm] = prms[obj][prm][0]
            
    prm_list, key_list = a.fit.flatten_prms(start_prms)
    prm_list = n.array(prm_list)
    #End: Parse Parameters
    



    #stream the data to the nodes
    print "streaming out data "
    print "fitting will occur as we work through the stream"
#    cat_head = {'srcs':['cyg'],'catalogs':['misc']}
#    comm.bcast(cat_head,root=0)
    tnow=0
    first=True
    data = n.array([],dtype=uv_aipy,ndmin=1)
    #print a nice header for the sample trace
    out = '#jd \t lst \t'
    for name,parms in prm_dict.iteritems(): 
        for parm in parms: out += name+'_'+parm+'\t'
        out += "\t"
    out += "score \t first_score\n"
    if not opts.quiet: print out
    if not opts.tracefile is None: 
        tracefile = open(opts.tracefile,'w')
        tracefile.write(out)
    
    for file in args:
        uv = a.miriad.UV(file)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        uv.select('decimate', opts.decimate, opts.decphs)
        for p,d in uv.all():
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
            if tnow !=t and first:
                tnow = t
                first = False
                rec = n.array((p,d,mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,rec))
            elif tnow==t:
                rec = n.array((p,d,mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,rec))                
            elif tnow != t and not first:
                #scatter the data, taking care not to send any to rank=0
                comm.bcast(2,root=0) #cmd everybody to load data
                comm.scatter([0]*1+n.array_split(data,size-1),root=0)
                print "cache summary"
                print '   %d samples'%n.sum(n.ones_like(data['spec']))
#                except TypeError: print '   %d samples'%len(data)
                print '   %d valid'%n.sum(n.logical_not(data['mask'].astype(n.int)))
                print '   %f sum weights'%n.sum(data['preamble']['m1'])
                #do the fit for this time
                first_fit=None
                rv = a.optimize.fmin(fit_func,prm_list,
                    full_output=1,disp=0,
                    maxfun=opts.maxiter,maxiter=n.Inf,
                    ftol=opts.ftol,xtol=opts.xtol)              
                if not opts.quiet: print str(t)+"\t"+print_trace(rv,prm_dict)
                prms,score = rv[:2]
                prms = a.fit.reconstruct_prms(prms, key_list)
                a.fit.print_params(prms)
                print 'Score:', score * first_fit, 
                print '(%2.2f%% of %f)' % (100 * float(score), first_fit)
                print '------------------------------------------------------------'
                if not opts.tracefile is None: 
                    tracefile.write(str(t)+'\t' +str(uv['lst']) + '\t'+print_trace(rv,prm_dict))
                    tracefile.flush()

                #get ready for the next time step
#                comm.bcast(0,root=0) #reset the other nodes
                tnow=t
                data = n.array([],dtype=uv_aipy,ndmin=1) #clear the time buffer
                rec = n.array((p,d,mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,rec))
    comm.bcast(0,root=0) #reset the other nodes
else: #All other nodes
#    cat_head = comm.bcast(None)
#    cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
    # sim difference loop
    print rank,"starting sim loop"
    data = None
    while(True):
        #decide what to do
        cmd = comm.bcast(root=0)
        if cmd==2: data=comm.scatter(root=0);continue #load data
        elif cmd==1: pass #continue on to fitting
        elif cmd==0: break #exit loop
        key_list = comm.bcast(None,root=0)
        prm_list = comm.bcast(None,root=0)
        #------
        
        prms = a.fit.reconstruct_prms(prm_list, key_list)
#        a.fit.print_params(prms)
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
                eqs = cat.get_crds('eq', ncrd=3)
                flx = cat.get_jys()
                dra,ddec = cat.get('ionref')
                aa.sim_cache(eqs, flx, mfreqs=mfq, 
                    ionrefs=(dra,ddec), srcshapes=(a1,a2,th))
            i,j = p['bl']
#                print '\t',i,j
            if i==j: continue
            simd = aa.sim(i,j,pol=p['pol'])
            difsq = n.abs(d - simd)**2
#            if i==10 and j==14:
#                print d,simd
#                print p['t'],i,j,difsq
#                comm.bcast(None,root=0)
            difsq = n.where(f,0,difsq)
            score += difsq.sum()
            m1 += p['m1']
        comm.gather(score)
        comm.gather(m1)