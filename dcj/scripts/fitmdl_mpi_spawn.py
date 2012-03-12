#!/usr/bin/env python
#
#  fitmdl_mpi_spawn.py
#  
#
#  Created by Danny Jacobs on 4/14/10.
#  PAPER Project
#
"""
A fitmdl called by a higher level fitter.
Options input as a list
args = [sdf,sfreq,nchan,opts.chan,opts.src,opts.cat,opts.cal]
Just enough to creat cat and aa objects, which for some dumb reason
aren't serializable.

Postpone until I can get analyzer results. This will only help if scatter is the 
bottleneck.
Wed 14 April 2010

"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
from mpi4py import MPI


comm = MPI.Comm.Get_parent()
size = comm.Get_size()
rank = comm.Get_rank()

sdf,sfreq,nchan,opts.chan,opts.src,opts.cat,opts.cal = 
nchan = uv['nchan']
aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
while(True)
    cmd = comm.bcast(root=0)
    if cmd==0:
        sys.exit()
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
        aa.set_params(prms)
        cat.set_params(prms)
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