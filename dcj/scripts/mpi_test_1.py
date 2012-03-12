#!/usr/bin/env python
#
#  mpi_test_1.py
#  
#
#  Created by Danny Jacobs on 3/23/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pickle,time
from mpi4py import MPI

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--test',type="int",
    help="Choose a test number.")
opts, args = o.parse_args(sys.argv[1:])


"""
Functions
"""
def about(uv):
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


"""
Test 1: numpy arrays
"""

if opts.test==1:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # automatic MPI datatype discovery
    # for NumPy arrays and PEP-3118 buffers
    if rank == 0:
       data = n.arange(100, dtype=n.float64)
       print "sending",n.max(data),data.shape
       comm.Send(data, dest=1, tag=13)
    elif rank == 1:
       data = n.empty(100, dtype=n.float64)
       comm.Recv(data, source=0, tag=13)
       print "recieved ", n.max(data), data.shape 
#    # pass explicit MPI datatypes
#    if rank == 0:
#       data = numpy.arange(1000, dtype='i')
#       comm.Send([data, MPI.INT], dest=1, tag=77)
#    elif rank == 1:
#       data = numpy.empty(1000, dtype='i')
#       comm.Recv([data, MPI.INT], source=0, tag=77)
if opts.test==2:
    "Broadcast test"
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank==0:
        data = {'list':[4,5,5],
                'tuple':('a','n')
            }
    else:
        data = None
    data = comm.bcast(data,root=0)
    print rank,data
if opts.test==3:
    "Broadcast/gather test"
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank==0:
        data = n.ones(10)
    else:
        data = None
    data = comm.bcast(data,root=0)
    data *= 2
    data = comm.allgather(data)
    if rank==0: print rank,data
if opts.test==4:
    "Reduce test"
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank==0:
        data = n.ones(1,'d')
    else:
        data = None    
    data = comm.bcast(data,root=0)
    data *=2
    result = comm.reduce(data,None,op=MPI.SUM,root=0)
    if rank==0: print rank, result,"is it 4? it should be 4"
if opts.test==5:
    "My pi computation: pympi example using mpi4py"
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    oldpi, pi, mypi = 0.0,0.0,0.0
    if rank==0:
        N = 100000
    else:
        N = None
    N = comm.bcast(N,root=0)
    inside=0
    for i in range(N):
        x = n.random.uniform(0,1)
        y = n.random.uniform(0,1)
        if ((x*x) + (y*y))<1:
            inside += 1
        oldpi = pi
        mypi = (inside*1.0)/N
        #print rank,mypi*4.0/size
        comm.reduce(mypi, None, op=MPI.SUM, root=0)
        pi = ((4.0) * mypi + pi)/2
    print "uberpi",n.average(comm.allgather(pi))
#    print "pi = %15.14f"%pi
if opts.test==6:
    "basic numpy array broadcast test"
    "basic buffer, sim, score, return"
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank==0:
        data = n.ones(10,dtype=n.float64)  #fake load the data
    else:
        data = n.zeros(10,dtype=n.float64)
    print rank,data
    comm.Bcast(data,root=0)
    print rank,data,'.'    
if opts.test==7:
    "basic buffer, sim, score, return"
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    ncount = 100
    if rank==0:
        dsize = 10
        data = n.ones(dsize,dtype=n.float64)  #fake load the data
        dsize = comm.bcast(dsize,root=0)
    else:
        dsize = None
        dsize = comm.bcast(dsize,root=0)
        data = n.zeros(dsize,dtype=n.float64)
    comm.Bcast(data,root=0)
    score =0
    for i in range(ncount):
        sim = n.random.normal(size=len(data))
        score += n.average(data-sim)
        comm.reduce(score,None,op=MPI.SUM,root=0)
    score /= ncount
    score= n.average(comm.allgather(score))
    if rank==0: print score
if opts.test==8:
    """
    Test "any object can be passed"
    ---> morphed into miriad file buffer!
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    uv = a.miriad.UV('test.uv')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2)]), 
                    ('spec', n.float,nchan),
                    ('mask',n.bool,nchan)])
    # automatic MPI datatype discovery
    # for NumPy arrays and PEP-3118 buffers
    if rank == 0:
        uv = a.miriad.UV('test.uv')
        data = n.array([],dtype=uv_aipy,ndmin=1)
        recnum = 0
        for p,d in uv.all():
            new = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
            data = n.concatenate((data,new))
        chunks = n.array_split(data,size-1)
        for i in range(1,size):
            comm.send(chunks[i-1], dest=i, tag=1)
    else:
        data = comm.recv(source=0,tag=1)
        print rank,"recieved ",len(data),"records"
if opts.test==9:
    "test sending metadata"
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    uv = a.miriad.UV('test.uv')
    head = {}
    for k in uv.vars():
        if k not in ['corr']: head[k] = uv[k]
    if rank==0:
        for i in range(1,size):
            comm.send(head,dest=i,tag=10)
    else:
        uv_head = comm.recv(source=0,tag=10)
        print uv_head
if opts.test==10:
    "can we scatter a list of arrays?"
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    uv = a.miriad.UV('new.uv')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2)]),
                    ('spec', n.complexfloating,nchan),
                    ('mask',n.bool,nchan)])
    if rank == 0:
        #buffer data completely before sending
        uv = a.miriad.UV('new.uv')
        data = n.array([],dtype=uv_aipy,ndmin=1)
        recnum = 0
        tnow =0
        first = True
        for p,d in uv.all():
            #get a single time

            if tnow != p[1] and first: 
                tnow=p[1]
                first = False
                new = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,new))
            elif tnow ==p[1]:
                new = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,new))            
            elif tnow != p[1] and not first:
                print "accumulated ",len(data),"records"
                print "scattering bls for time =",p[1]
                print "scattering them!"
                comm.scatter(n.concatenate([[0]*1,n.array_split(data,size-1)]),root=0)
                break

            
    else:
        data = comm.scatter(0,root=0)
        print rank,"recieved ",len(data),"records"
        print rank,"recieved  baslines:",data['preamble']['bl']
    #yes we can! (but we can't scatter just any iterable eg uv.all())
if opts.test==11:
    "test stream scattering data"
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    uv = a.miriad.UV('test.uv')
    def about(uv):
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
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2)]), 
                    ('spec', n.float,nchan),
                    ('mask',n.bool,nchan)])
    if rank==0:
        #divide along t only
        tstart = time.time()
        ts,bls,freqs,nrec = about(uv)
        print rank, "read",nrec,"records"
        print rank, "with ",len(ts),"times and ",len(bls),"baselines"
        tend = time.time()
        print "times",tstart,tend
        print "IO speed", nrec/(tend-tstart)," [recs/second]"
        tset = n.array_split(ts,size-1)
        nsent =n.zeros(size-1)
        for p,d in uv.all():
            rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
            #find the dest node
            for i in range(1,size):
                if rec['preamble']['t'] in tset[i-1]:
                    comm.bcast(i,root=0)
                    comm.send(rec, dest=i, tag=1)
                    nsent[i-1] += 1
        comm.bcast(0,root=0)
        print rank, "node allocations",nsent
        print rank, "total sent", n.sum(nsent)
    else:
        data = n.array([],dtype=uv_aipy,ndmin=1)
        while(True):
            inc_rank = comm.bcast(root=0)
            if inc_rank==0: print "breaking";break
            elif inc_rank==rank:
                _in = comm.recv(source=0,tag=1)
                data = n.concatenate((data,_in))
#                print rank,len(data),
        print rank,"recieved ",len(data),"records"
if opts.test==12:
    "test stream scattering data [condensed version of above]"
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    uv = a.miriad.UV('test.uv')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2)]), 
                    ('spec', n.float,nchan),
                    ('mask',n.bool,nchan)])
    if rank==0:  #The "master node"
        #divide along t only
        ts,bls,freqs,nrec = about(uv)
        tset = n.array_split(ts,size-1)
        #stream the data to the nodes
        for p,d in uv.all():
            rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
            for i in range(1,size):
                if rec['preamble']['t'] in tset[i-1]:
                    comm.bcast(i,root=0)               #alert the grid about incoming
                    comm.send(rec, dest=i, tag=1)      #send the record
        comm.bcast(0,root=0)                           #alert the grid: finished buffering
    else: #All other nodes
        data = n.array([],dtype=uv_aipy,ndmin=1)      #Create the data container
        while(True):
            inc_rank = comm.bcast(root=0)             #Recieve instructions
            if inc_rank==0: break                     #Quit if we're done
            elif inc_rank==rank:                      #otherwise, recieve data, if its for us
                _in = comm.recv(source=0,tag=1)
                data = n.concatenate((data,_in))
        print rank,"recieved ",len(data),"records"
if opts.test==13:
    "Test sim setup"
    #buffer data (do I _have_ to do this?)
    #send array, cal, and cat
    #compute sim
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    uv = a.miriad.UV('test.uv')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2)]), 
                    ('spec', n.float,nchan),
                    ('mask',n.bool,nchan)])
    if rank==0:
        uv = a.miriad.UV('test.uv')
        head = {}
        for k in uv.vars():
            if k not in ['corr']: head[k] = uv[k]
        head['cal']='pgb015_v005_gc'
        comm.send(head,dest=1,tag=1)
        cat_head = {'srcs':['286'],'catalogs':['three_cr']}
        comm.send(cat_head,dest=1,tag=2)
        t = 2455015.5
#        aa.set_jultime(t)
        comm.send(t,dest=1,tag=3)
    elif rank==1:
        uv_head = comm.recv(source=0,tag=1)
        aa = a.cal.get_aa(uv_head['cal'],uv_head['sdf'],uv_head['sfreq'],uv_head['nchan'])
        cat_head = comm.recv(source=0,tag=2)
        cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
        t = comm.recv(source=0,tag=3)
        aa.set_jultime(t)  
        cat.compute(aa)
        print len(aa),len(cat.get_jys()),t
        aa.sim_cache(cat.get_crds('eq',ncrd=3),cat.get_jys(),
            mfreqs=cat.get('mfreq'))
        print rank, "simulating..."
        simd = aa.sim(0,1,pol='yy')
        print rank,len(simd),n.max(simd),n.min(simd)
    else:
        print "Nothing!"
if opts.test==14:
    "basic fitmdl (no options)"
    "features: stream scattering, automatic time division"
    "no catalog division, only works with n=2"
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    uv = a.miriad.UV('test.uv')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2),
                         ('m1',n.complexfloating)]),
                    ('spec', n.float,nchan),
                    ('mask',n.bool,nchan)])
    if rank==0:  #The "master node"
        #divide along t only
        ts,bls,freqs,nrec = about(uv)
        tset = n.array_split(ts,size-1)
        #stream the data to the nodes
        for p,d in uv.all():
            m2 = n.where(d.mask,0,n.ma.abs(d)**2).sum()
            p += (m2,)
            rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
            for i in range(1,size):
                if rec['preamble']['t'] in tset[i-1]:
                    comm.bcast(i,root=0)               #alert the grid about incoming
                    comm.send(rec, dest=i, tag=1)      #send the record
        comm.bcast(0,root=0)                           #alert the grid: finished buffering
        #set up the AntennaArray and Catalog
        uv.rewind()
        uv_head = {}
        for k in uv.vars():
            if k not in ['corr']: uv_head[k] = uv[k]
        uv_head['cal']='pgb015_v005_gc'
        comm.send(uv_head,dest=1,tag=13)
        cat_head = {'srcs':['286'],'catalogs':['three_cr']}
        comm.send(cat_head,dest=1,tag=14)
        score = 0
        m1 = 0
        scores = comm.allgather(score)
        m1s = comm.allgather(m1)
        print "number of scores = %d, score = %f, m1 = %f"%(len(scores),n.sum(scores)/n.sum(m1s),n.sum(m1s))
    else: #All other nodes
        data = n.array([],dtype=uv_aipy,ndmin=1)      #Create the data container
        #Buffer loop
        while(True):
            inc_rank = comm.bcast(root=0)             #Recieve instructions
            if inc_rank==0: break                     #Quit if we're done
            elif inc_rank==rank:                      #otherwise, recieve data, if its for us
                _in = comm.recv(source=0,tag=1)
                data = n.concatenate((data,_in))
        print rank,"recieved ",len(data),"records"
        #set up the AntennaArray,Catalog
        uv_head = comm.recv(source=0,tag=13)
        aa = a.cal.get_aa(uv_head['cal'],uv_head['sdf'],uv_head['sfreq'],uv_head['nchan'])
        cat_head = comm.recv(source=0,tag=14)
        cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
        # sim difference loop
        print "starting sim loop"
        score = 0
        m1 = 0
        for p,d,f in data:
            aa.set_jultime(p['t'])
            cat.compute(aa)
            aa.sim_cache(cat.get_crds('eq',ncrd=3),cat.get_jys(),
                mfreqs=cat.get('mfreq'))
            simd = aa.sim(0,1,pol='yy')
            difsq = n.abs(d - simd)**2
            difsq = n.where(f,0,difsq)
            score += difsq.sum()
            m1 += p['m1']
        comm.allgather(score)
        comm.allgather(m1)
if opts.test==15:
    "basic score calc (no options)"
    "features: stream scattering"
    "send the array and catalog parameters to everyone"
    "get a score back"
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
#    uv = a.miriad.UV('test.uv')
    print "opening ",'zen.2454565.41207.uvc'
    uv = a.miriad.UV('zen.2454565.41207.uvc')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2),
                         ('m1',n.complexfloating),
                         ('pol',str,2)]),
                    ('spec', n.float,nchan),
                    ('mask',n.bool,nchan)])
    if rank==0:  #The "master node"
        #divide along t only
        ts,bls,freqs,nrec = about(uv)
        tset = n.array_split(ts,size-1)
        #stream the data to the nodes
        print "streaming out data"
        for p,d in uv.all():
            print '.',
            m2 = n.where(d.mask,0,n.ma.abs(d)**2).sum()
            p += (m2,)
            p += (a.miriad.pol2str[uv['pol']],)
            rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
            for i in range(1,size):
                if rec['preamble']['t'] in tset[i-1]:
                    comm.bcast(i,root=0)               #alert the grid about incoming
                    comm.send(rec, dest=i, tag=1)      #send the record
        comm.bcast(0,root=0)                           #alert the grid: finished buffering
        #set up the AntennaArray and Catalog
        uv.rewind()
        uv_head = {}
        for k in uv.vars():
            if k not in ['corr']: uv_head[k] = uv[k]
#        uv_head['cal']='pgb015_v005_gc'
        uv_head['cal'] = 'pgb564'
        comm.bcast(uv_head,root=0)
        cat_head = {'srcs':['286'],'catalogs':['three_cr']}
        comm.bcast(cat_head,root=0)
        score = 0
        m1 = 0
        scores = comm.allgather(score)
        m1s = comm.allgather(m1)
        print "number of scores = %d, score = %f, m1 = %f"%(len(scores),n.sum(scores)/n.sum(m1s),n.sum(m1s))
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
        print rank,"recieved ",len(data),"records"
        #set up the AntennaArray,Catalog
        uv_head = comm.bcast(None)
        aa = a.cal.get_aa(uv_head['cal'],uv_head['sdf'],uv_head['sfreq'],uv_head['nchan'])
        cat_head = comm.bcast(None)
        cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
        # sim difference loop
        print "starting sim loop"
        score = 0
        m1 = 0
        for p,d,f in data:
            aa.set_jultime(p['t'])
            cat.compute(aa)
            aa.sim_cache(cat.get_crds('eq',ncrd=3),cat.get_jys(),
                mfreqs=cat.get('mfreq'))
            simd = aa.sim(0,1,pol=p['pol'])
            difsq = n.abs(d - simd)**2
            difsq = n.where(f,0,difsq)
            score += difsq.sum()
            m1 += p['m1']
        comm.allgather(score)
        comm.allgather(m1)
if opts.test==16:
    "basic fitmdl (no options)"
    "features: stream scattering"
    "send the array and catalog parameters to everyone"
    "send new parameters a few times in a loop"
    #Begin: Setup MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    #End: Setup MPI 

    #Begin: Setup environment
    uv = a.miriad.UV('new.uv')
    incal = 'pgb015_v005_gc'
    inshprms = None
    inprms = "cyg=jys"
    srcs = ['cyg']
    cats = ['misc']
#    print "opening ",'zen.2454565.41207.uvc'
#    uv = a.miriad.UV('zen.2454565.41207.uvc')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2),
                         ('m1',n.complexfloating),
                         ('pol',str,2)]),
                    ('spec', n.float,nchan),
                    ('mask',n.bool,nchan)])
    aa = a.cal.get_aa(incal,uv['sdf'],uv['sfreq'],uv['nchan'])
    cat = a.src.get_catalog(srcs=srcs,catalogs=cats)
    if rank==0:  #The "master node"
        #Begin: Parse Parameters
        prms, prm_dict, shkeys = {}, {}, []
        # Handle shared parameters
        if inshprms:
            shprms = map(a.scripting.parse_prms, inshprms.split(','))
            for s in shprms:
                keys = s.keys(); keys.sort()
                k = keys[0]
                # Only enter one instance of shared parameters (indexed under first key)
                if prms.has_key(k): prms[k].update(s[k])
                else: prms.update({k:s[k]})
                # Add an entry to shkeys for propagating variables to other objects
                shkeys.append((keys, s[k].keys()))
        # Handle normal parameters
        if not inprms is None:
            pd = a.scripting.parse_prms(inprms)
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
        


    #End: Setup environment
        #divide along t only
        ts,bls,freqs,nrec = about(uv)
        tset = n.array_split(ts,size-1)
        #stream the data to the nodes
        print "streaming out data"
        for p,d in uv.all():
            m2 = n.where(d.mask,0,n.ma.abs(d)**2).sum()
            p += (m2,)
            p += (a.miriad.pol2str[uv['pol']],)
            rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
            for i in range(1,size):
                if rec['preamble']['t'] in tset[i-1]:
                    comm.bcast(i,root=0)               #alert the grid about incoming
                    comm.send(rec, dest=i, tag=1)      #send the record
        comm.bcast(0,root=0)                           #alert the grid: finished buffering
        #set up the AntennaArray and Catalog
        uv.rewind()
        cat_head = {'srcs':['cyg'],'catalogs':['misc']}
        comm.bcast(cat_head,root=0)
        score = 0
        m1 = 0
        for j in range(2):
            comm.bcast(1,root=0)
            comm.bcast(key_list,root=0)
            comm.bcast(prm_list,root=0)
            scores = comm.gather(score)
            m1s = comm.gather(m1)
            print "number of scores = %d, score = %f, m1 = %f"%(len(scores),n.sum(scores)/n.sum(m1s),n.sum(m1s))
        comm.bcast(0,root=0)
    else: #All other nodes
        data = n.array([],dtype=uv_aipy,ndmin=1)      #Create the data container
        #Buffer loop
        while(True):
            inc_rank = comm.bcast(root=0)             #Recieve instructions
            if inc_rank==0: break                     #Quit if we're done
            elif inc_rank==rank:                      #otherwise, recieve data, if its for us
                _in = comm.recv(source=0,tag=1)
                data = n.concatenate((data,_in))
        print rank,"recieved ",len(data),"records"
        cat_head = comm.bcast(None)
        cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
        # sim difference loop
        print rank,"starting sim loop"
        while(True):
            if comm.bcast(root=0)==0: break
            key_list = comm.bcast(None,root=0)
            prm_list = comm.bcast(None,root=0)
            prms = a.fit.reconstruct_prms(prm_list, key_list)
            #a.fit.print_params(prms)
            aa.set_params(prms)
            cat.set_params(prms)
            score = 0
            m1 = 0
            for p,d,f in data:
                aa.set_jultime(p['t'])
                cat.compute(aa)
                aa.sim_cache(cat.get_crds('eq',ncrd=3),cat.get_jys(),
                    mfreqs=cat.get('mfreq'))
                simd = aa.sim(0,1,pol=p['pol'])
                difsq = n.abs(d - simd)**2
                difsq = n.where(f,0,difsq)
                score += difsq.sum()
                m1 += p['m1']
            comm.gather(score)
            comm.gather(m1)
if opts.test==17:
    "basic fitmdl (no options)"
    "features: stream scattering"
    "send the array and catalog parameters to everyone"
    "fit averaging overall deps (t,bl,freqs)"
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
        if first_fit is None: 
            first_fit = score
            return score
        print "number of scores = %d, score = %f, m1 = %f"%(len(scores),n.sum(scores)/n.sum(m1s),n.sum(m1s))
        return score/first_fit

    
    #Begin: Setup environment
    uv = a.miriad.UV('new.uv')
    a.scripting.uv_selector(uv, 'cross', 'xx')
    incal = 'pgb015_v005'
    inshprms = None
    inprms = "cyg=jys"
    srcs = ['cyg']
    cats = ['misc']
#    print "opening ",'zen.2454565.41207.uvc'
#    uv = a.miriad.UV('zen.2454565.41207.uvc')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2),
                         ('m1',n.complexfloating),
                         ('pol',str,2)]),
                    ('spec', n.complexfloating,nchan),
                    ('mask',n.bool,nchan)])
    aa = a.cal.get_aa(incal,uv['sdf'],uv['sfreq'],uv['nchan'])
    cat = a.src.get_catalog(srcs=srcs,catalogs=cats)
    if rank==0:  #The "master node"
        #Begin: Parse Parameters
        prms, prm_dict, shkeys = {}, {}, []
        # Handle shared parameters
        if inshprms:
            shprms = map(a.scripting.parse_prms, inshprms.split(','))
            for s in shprms:
                keys = s.keys(); keys.sort()
                k = keys[0]
                # Only enter one instance of shared parameters (indexed under first key)
                if prms.has_key(k): prms[k].update(s[k])
                else: prms.update({k:s[k]})
                # Add an entry to shkeys for propagating variables to other objects
                shkeys.append((keys, s[k].keys()))
        # Handle normal parameters
        if not inprms is None:
            pd = a.scripting.parse_prms(inprms)
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
        


    #End: Setup environment
        #divide along t only
        ts,bls,freqs,nrec = about(uv)
        tset = n.array_split(ts,size-1)
        #stream the data to the nodes
        print "streaming out data "
        for p,d in uv.all():
            m2 = n.where(d.mask,0,n.ma.abs(d)**2).sum()
            p += (m2,)
            p += (a.miriad.pol2str[uv['pol']],)
            rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
            for i in range(1,size):
                if rec['preamble']['t'] in tset[i-1]:
                    comm.bcast(i,root=0)               #alert the grid about incoming
                    comm.send(rec, dest=i, tag=1)      #send the record
        comm.bcast(0,root=0)                           #alert the grid: finished buffering
        #set up the AntennaArray and Catalog
        uv.rewind()
        cat_head = {'srcs':['cyg'],'catalogs':['misc']}
        comm.bcast(cat_head,root=0)
        rv = a.optimize.fmin(fit_func,prm_list,
            full_output=1,disp=0,
            maxfun=n.Inf,maxiter=n.Inf,
            ftol=1e-10,xtol=1e-10)
        comm.bcast(0,root=0)#clear the nodes
    else: #All other nodes
        data = n.array([],dtype=uv_aipy,ndmin=1)      #Create the data container
        #Buffer loop
        while(True):
            inc_rank = comm.bcast(root=0)             #Recieve instructions
            if inc_rank==0: break                     #Quit if we're done
            elif inc_rank==rank:                      #otherwise, recieve data, if its for us
                _in = comm.recv(source=0,tag=1)
                data = n.concatenate((data,_in))
        print rank,"recieved ",len(data),"records having",len(data)*len(data[0]['spec']),"samples"
        cat_head = comm.bcast(None)
        cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
        # sim difference loop
        print rank,"starting sim loop"
        while(True):
            if comm.bcast(root=0)==0: break
            key_list = comm.bcast(None,root=0)
            prm_list = comm.bcast(None,root=0)
            prms = a.fit.reconstruct_prms(prm_list, key_list)
            a.fit.print_params(prms)
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
if opts.test==18:
    "basic fitmdl (no options)"
    "features: stream scattering"
    "send the array and catalog parameters to everyone"
    "fit averaging over deps (bl,freqs), leaving t independant"
    "'snap' mode, a solution for every time"
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
        if first_fit is None: 
            first_fit = score
            return score
        print "number of scores = %d, score = %f, m1 = %f"%(len(scores),n.sum(scores)/n.sum(m1s),n.sum(m1s))
        return score/first_fit

    
    #Begin: Setup environment
    uv = a.miriad.UV('new.uv')
    a.scripting.uv_selector(uv, 'cross', 'xx')
    incal = 'pgb015_v005'
    inshprms = None
    inprms = "cyg=jys"
    srcs = ['cyg']
    cats = ['misc']
#    print "opening ",'zen.2454565.41207.uvc'
#    uv = a.miriad.UV('zen.2454565.41207.uvc')
    nchan = uv['nchan']
    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                         ('t', n.float), 
                         ('bl', '<i4', 2),
                         ('m1',n.complexfloating),
                         ('pol',str,2)]),
                    ('spec', n.complexfloating,nchan),
                    ('mask',n.bool,nchan)])
    aa = a.cal.get_aa(incal,uv['sdf'],uv['sfreq'],uv['nchan'])
    cat = a.src.get_catalog(srcs=srcs,catalogs=cats)
    if rank==0:  #The "master node"
        #Begin: Parse Parameters
        prms, prm_dict, shkeys = {}, {}, []
        # Handle shared parameters
        if inshprms:
            shprms = map(a.scripting.parse_prms, inshprms.split(','))
            for s in shprms:
                keys = s.keys(); keys.sort()
                k = keys[0]
                # Only enter one instance of shared parameters (indexed under first key)
                if prms.has_key(k): prms[k].update(s[k])
                else: prms.update({k:s[k]})
                # Add an entry to shkeys for propagating variables to other objects
                shkeys.append((keys, s[k].keys()))
        # Handle normal parameters
        if not inprms is None:
            pd = a.scripting.parse_prms(inprms)
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
        


    #End: Setup environment
        #divide along t only
        ts,bls,freqs,nrec = about(uv)
        tset = n.array_split(ts,size-1)
        #stream the data to the nodes
        print "streaming out data "
        print "fitting will occur as we work through the stream"
        cat_head = {'srcs':['cyg'],'catalogs':['misc']}
        comm.bcast(cat_head,root=0)
        tnow=0
        first=True
        data = n.array([],dtype=uv_aipy,ndmin=1)
        for p,d in uv.all():
            m2 = n.where(d.mask,0,n.ma.abs(d)**2).sum()
            p += (m2,)
            p += (a.miriad.pol2str[uv['pol']],)
            t = p[1]
            if tnow !=t and first:
                tnow = t
                first = False
                rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,rec))
            elif tnow==t:
                rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,rec))                
            elif tnow != t and not first:
                #scatter the data, taking care not to send any to rank=0
                comm.bcast(2,root=0) #cmd everybody to load data
                comm.scatter([0]*1+n.array_split(data,size-1),root=0)
                #do the fit for this time
                rv = a.optimize.fmin(fit_func,prm_list,
                    full_output=1,disp=0,
                    maxfun=n.Inf,maxiter=n.Inf,
                    ftol=1e-10,xtol=1e-10)                
                #get ready for the next time step
                comm.bcast(0,root=0) #reset the other nodes
                tnow=t
                data = n.array([],dtype=uv_aipy,ndmin=1) #clear the time buffer
                rec = n.array((p,d.data,d.mask),dtype=uv_aipy,ndmin=1)
                data = n.concatenate((data,rec))
    else: #All other nodes
        cat_head = comm.bcast(None)
        cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
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
            a.fit.print_params(prms)
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
if opts.test==19:
    "Optimize variable catalogs. Can I really truly not serialize them?"
    #buffer data (do I _have_ to do this?)
    #send array, cal, and cat
    #compute sim
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
#    uv = a.miriad.UV('test.uv')
#    nchan = uv['nchan']
#    uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
#                         ('t', n.float), 
#                         ('bl', '<i4', 2)]), 
#                    ('spec', n.float,nchan),
#                    ('mask',n.bool,nchan)])
    print "node",rank,"starting"
    if rank==0:
#        uv = a.miriad.UV('test.uv')
#        head = {}
#        for k in uv.vars():
#            if k not in ['corr']: head[k] = uv[k]
#        head['cal']='pgb015_v005_gc'
#        comm.send(head,dest=1,tag=1)
#        cat_head = {'srcs':['286'],'catalogs':['three_cr']}
        fullcat = a.src.get_catalog(catalogs=['three_cr'])
        for src in fullcat.keys():
            if type(fullcat[src])==a.fit.RadioSpecial:
                fullcat.pop(src)
        print fullcat.keys()
        #cat = a.src.get_catalog(catalogs=['misc'],srcs=srcs[1:3])
        print "node 0 catalog loaded, Sending..."
        comm.send(fullcat,dest=1,tag=1)
        aa = a.cal.get_aa('pgb015_v005',0.00005,0.13,10)
        comm.send(aa,dest=1,tag=2)
        
#        t = 2455015.5
##        aa.set_jultime(t)
#        comm.send(t,dest=1,tag=3)
    elif rank==1:
        print "waiting for catalog from node 0"
        cat = comm.recv(source=0,tag=1)
        print cat
        aa = comm.recv(source=0,tag=2)
        
#        uv_head = comm.recv(source=0,tag=1)
#        aa = a.cal.get_aa(uv_head['cal'],uv_head['sdf'],uv_head['sfreq'],uv_head['nchan'])
#        cat_head = comm.recv(source=0,tag=2)
#        cat = a.src.get_catalog(srcs=cat_head['srcs'],catalogs=cat_head['catalogs'])
#        t = comm.recv(source=0,tag=3)
#        aa.set_jultime(t)  
#        cat.compute(aa)
#        print len(aa),len(cat.get_jys()),t
#        aa.sim_cache(cat.get_crds('eq',ncrd=3),cat.get_jys(),
#            mfreqs=cat.get('mfreq'))
#        print rank, "simulating..."
#        simd = aa.sim(0,1,pol='yy')
#        print rank,len(simd),n.max(simd),n.min(simd)
    else:
        print "Nothing!"        