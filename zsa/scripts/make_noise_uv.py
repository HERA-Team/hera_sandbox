#! /usr/bin/env python
import aipy as a
import numpy as n
import optparse, sys
import capo 

aa = a.cal.get_aa('psa6240_v003', n.linspace(.1,.2,203))
seps, conj = capo.red.group_redundant_bls(aa.ant_layout)

o = optparse.OptionParser()
o.add_option('--signal', action='store_true',
            help='Only store signal.')
o.add_option('--noise', action='store_true',
            help='Only store noise.')
opts,args = o.parse_args(sys.argv[1:])

curtime=None
noise_t_fr=None

def noise(nsamp, var):
    return n.random.normal(size=nsamp)*n.sqrt(var/2.) + 1j*n.random.normal(size=nsamp)*n.sqrt(var/2. )

def mfunc(uv, p, d, f):
    global curtime, noise_t_fr
    uvw, t, (i,j) = p
    bl = a.miriad.ij2bl(i,j)
    if curtime != t:
        curtime = t
        #new noise for every time step/freq but same for all baselines.
        #noise_t_fr = n.random.normal(size=len(d)) * n.exp(2j*n.pi*n.random.uniform(size=len(d)))
        noise_t_fr = noise(len(d), 1.0)
    #noise_t_fr_bl = n.random.normal(size=len(d)) * n.exp(2j*n.pi*n.random.uniform(size=len(d)))
    noise_t_fr_bl = noise(len(d), 1.0)

    if opts.signal:
        if conj[bl]:
            data = n.conj(noise_t_fr) 
        else:
            data = noise_t_fr
        
    elif opts.noise:
        if conj[bl]:
            data = n.conj(noise_t_fr_bl)
        else:
            data = noise_t_fr_bl
    
    else:
        if conj[bl]:
            data = n.conj(noise_t_fr) + n.conj(noise_t_fr_bl)
        else:
            data = noise_t_fr + noise_t_fr_bl
    
    return p, data, f

for filename in args:
    filename_end = filename.split('/')[-1]
    if opts.signal:
        outfile = filename_end + '_signal'
    elif opts.noise:
        outfile = filename_end + '_noise'
    else:
        outfile = filename_end + 's+n'
     
    print 'Writing %s'%(outfile)
    uvi = a.miriad.UV(filename)
#    a.scripting.uv_selector(uvi, ants=baselines)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+'\n', raw=True)
