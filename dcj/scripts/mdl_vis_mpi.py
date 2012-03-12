#!/usr/bin/env python
#
#  mdl_dev_1.py
#  
#
#  Created by Danny Jacobs on 6/25/10.
#  PAPER Project
#

import aipy as a, numpy as n, math as m,os
import sys, optparse,pickle,time,curses as C,logging,warnings
from mpi4py import MPI
import tabular as tab,pickle


"""
A 40% new mpi-based implimentation of mdlvis.

Currently only supports generation of new simulations.
"""

o = optparse.OptionParser()
o.set_usage('mdl_vis_mpi.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('-m','--mode', dest='mode', default='sim',
    help='Operation mode.  Can be "sim" (output simulated data), "sub" (subtract from input data), or "add" (add to input data).  Default is "sim"')
o.add_option('-f', '--flag', dest='flag', action='store_true',
    help='If outputting a simulated data set, mimic the data flagging of the original dataset.')
#o.add_option('-m', '--map', dest='map',
#    help='The Healpix map to use for simulation input.')
#o.add_option('--iepoch', dest='iepoch', default=ephem.J2000, 
#    help='The epoch of coordinates in the map. Default J2000.')
#o.add_option('--freq', dest='freq', default=.150, type='float',
#    help='Frequency of flux data in map.')
o.add_option('-n', '--noiselev', dest='noiselev', default=0., type='float',
    help='RMS amplitude of noise (Jy) added to each UV sample of simulation.')
o.add_option('--nchan', dest='nchan', default=256, type='int',
    help='Number of channels in simulated data if no input data to mimic.  Default is 256')
o.add_option('--sfreq', dest='sfreq', default=.075, type='float',
    help='Start frequency (GHz) in simulated data if no input data to mimic.  Default is 0.075')
o.add_option('--sdf', dest='sdf', default=.150/256, type='float',
    help='Channel spacing (GHz) in simulated data if no input data to mimic.  Default is .150/256')
o.add_option('--inttime', dest='inttime', default=10, type='float',
    help='Integration time (s) in simulated data if no input data to mimic.  Default is 10')
o.add_option('--startjd', dest='startjd', default=2454600., type='float',
    help='Julian Date to start observation if no input data to mimic.  Default is 2454600')
o.add_option('--endjd', dest='endjd', default=2454601., type='float',
    help='Julian Date to end observation if no input data to mimic.  Default is 2454601')
o.add_option('--pol', dest='pol', 
    help='Polarizations to simulate (xx,yy,xy,yx) if starting file from scratch.')
o.add_option('--prnt',
    help='comma delimited list of parameters to print. Only works with -m sim. NB for now only prints the header, you still have to make your cal file do the printing. ')
o.add_option('--sim_flags',type='float',
    help='Flag this fraction of a spectrum randomly to simulate rfi excision')
o.add_option('--out',default='new.uv',
    help='Output file name for new sims. [new.uv]')
opts, args = o.parse_args(sys.argv[1:])
logging.basicConfig(level=logging.ERROR)
log = logging.getLogger('mdl_dev_1')

def init_uv(uvofile):
    uv = a.miriad.UV(uvofile, status='new')
    uv._wrhd('obstype','mixed-auto-cross')
    uv._wrhd('history','MDLVIS_MPI: created file.\nMDLVIS: srcs=%s cat=%s mode=%s flag=%s noise=%f\n' % (opts.src, opts.cat, opts.mode, opts.flag, opts.noiselev))
    uv.add_var('telescop','a'); uv['telescop'] = 'PAPER'
    uv.add_var('operator','a'); uv['operator'] = 'AIPY'
    uv.add_var('version' ,'a'); uv['version'] = '0.0.1'
    uv.add_var('epoch'   ,'r'); uv['epoch'] = 2000.
    uv.add_var('source'  ,'a'); uv['source'] = 'zenith'
    uv.add_var('latitud' ,'d'); uv['latitud'] = aa.lat
    uv.add_var('dec'     ,'d'); uv['dec'] = aa.lat
    uv.add_var('obsdec'  ,'d'); uv['obsdec'] = aa.lat
    uv.add_var('longitu' ,'d'); uv['longitu'] = aa.long
    uv.add_var('npol'    ,'i'); uv['npol'] = len(pols)
    uv.add_var('nspect'  ,'i'); uv['nspect'] = 1
    uv.add_var('nants'   ,'i'); uv['nants'] = len(aa)
    uv.add_var('antpos'  ,'d')
    antpos = n.array([ant.pos for ant in aa], dtype=n.double)
    uv['antpos'] = antpos.transpose().flatten()
    uv.add_var('sfreq'   ,'d'); uv['sfreq'] = opts.sfreq
    uv.add_var('freq'    ,'d'); uv['freq'] = opts.sfreq
    uv.add_var('restfreq','d'); uv['freq'] = opts.sfreq
    uv.add_var('sdf'     ,'d'); uv['sdf'] = opts.sdf
    uv.add_var('nchan'   ,'i'); uv['nchan'] = opts.nchan
    uv.add_var('nschan'  ,'i'); uv['nschan'] = opts.nchan
    uv.add_var('inttime' ,'r'); uv['inttime'] = float(opts.inttime)
    # These variables just set to dummy values
    uv.add_var('vsource' ,'r'); uv['vsource'] = 0.
    uv.add_var('ischan'  ,'i'); uv['ischan'] = 1
    uv.add_var('tscale'  ,'r'); uv['tscale'] = 0.
    uv.add_var('veldop'  ,'r'); uv['veldop'] = 0.
    # These variables will get updated every spectrum
    uv.add_var('coord'   ,'d')
    uv.add_var('time'    ,'d')
    uv.add_var('lst'     ,'d')
    uv.add_var('ra'      ,'d')
    uv.add_var('obsra'   ,'d')
    uv.add_var('baseline','r')
    uv.add_var('pol'     ,'i')
    return uv

#Begin: Setup MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#End: Setup MPI
print rank,"Starting up"
def tracelog(file,prms,jd,lst):
    pflat,kys = a.fit.flatten_prms(prms)
    if not os.path.exists(file+'mdl_keys.pkl'):
        pickle.dump(kys,open(file+'mdl_keys.pkl','w'))
    open(file+'mdl_prms.txt','a').write('\t'.join(map(str,pflat))+'\n')
    open(file+'mdl_times.txt','a').write(str(jd)+'\t'+str(repr(lst))+'\n')

assert(len(args) > 0 or (opts.mode == 'sim' and not opts.flag and not (opts.pol is None)))
# Parse command-line options
if len(args) > 0:
    print "sorry, input files are not yet supported"
    sys.exit()
#    uv = a.miriad.UV(args[0])
#    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
#    p,d,f = uv.read(raw=True)
#    no_flags = n.zeros_like(f)
#    del(uv)
else:
    aa = a.cal.get_aa(opts.cal, opts.sdf, opts.sfreq, opts.nchan)
    no_data = n.zeros(opts.nchan, dtype=n.complex64)
    no_flags = n.zeros(opts.nchan, dtype=n.int32)
pols = opts.pol.split(',')

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
mfq = cat.get('mfreq')
a1s,a2s,ths = cat.get('srcshape')

#setup the data types for mpi data exchange
uv_aipy = n.dtype([('preamble', [('uvw', n.float, 3),
                     ('t', n.float), 
                     ('bl', '<i4', 2),
                     ('pol',str,2),
                     ('lst', n.float)]),
                ('spec', n.complexfloating,len(no_data)),
                ('mask',n.bool,len(no_data))])

# generate a vector of srcs, times and baselines
times = n.arange(opts.startjd, opts.endjd, opts.inttime/a.const.s_per_day)
bls = []
for i in range(len(aa)):
    for j in range(i,len(aa)):
        if i==j:continue
        bls.append([i,j])
bl_times = []
for t in times:
    for bl in bls:
        bl_times.append((bl[0],bl[1],t))
bl_times = comm.scatter(n.array_split(bl_times,size))
curtime = 0
data = n.array([],dtype=uv_aipy,ndmin=1)
datas = n.array([],dtype=uv_aipy,ndmin=1)
uvofile=opts.out
if os.path.exists(uvofile): 
    if rank==0: print uvofile,"exists, exiting..."
    sys.exit()
if rank != 0:
    for N,blt in enumerate(bl_times):
        (i,j),t = (map(int,blt[:2]),blt[-1])
        if i==j: continue
        if curtime!=t:
            curtime = t
    #        if opts.prnt is None: print rank,t,data.size
            aa.set_jultime(t)
            lst = aa.sidereal_time()
            cat.compute(aa)
     #       tracelog(uvofile+'/',cat.get_params().update(aa.get_params()),t,aa.sidereal_time()) 
            eqs = cat.get_crds('eq', ncrd=3)
            flx = cat.get_jys()
            sys.stdout.flush()
            dras, ddecs = cat.get('ionref')
            aa.sim_cache(eqs, flx, mfreqs=mfq, 
                ionrefs=(dras,ddecs), srcshapes=(a1s,a2s,ths))
        ai,aj = aa[i],aa[j]
        crd = ai.pos - aj.pos
        preamble = (crd, t, (i,j))
        for pol in pols:
            d = aa.sim(i, j, pol=pol)
            sys.stdout.flush()
            f = n.zeros_like(d).astype(n.int)
            preamble += (pol,)
            preamble += (lst,)
            rec = n.array((preamble,d,f),dtype=uv_aipy,ndmin=1)
            try: data = n.concatenate((data,rec))
            except(MemoryError): print rank,"MemoryError at concatenate";
        if data.size>1e2:
#            print rank,"Data size limit reached, gathering and concatenating"
            sys.stdout.flush()
            comm.gather(data.size)
            comm.gather(data,root=0)
            comm.gather(1)
            data = n.array([],dtype=uv_aipy,ndmin=1)
    #send out the remainder of the data
    comm.gather(data.size,root=0)
    comm.gather(data,root=0)
    comm.gather(0)
    print rank, "exiting"
else:
    cnt=0
    maxcnt=7e4
    exited = 0
    while True:
        print "gathering ",n.sum(comm.gather(0))," numbers"
        datas = comm.gather([])
        cont = comm.gather(1)
        print " %i nodes are on their last round"%(size-int(n.sum(cont)))
        data = n.concatenate(datas[1:])
        for p,d,f in data:
            if cnt%maxcnt==0:
                if cnt>0:del(uv)
                filename = '.'.join(('sim',str(p['t']),'uv'))
                print "creating: "+filename
                if os.path.exists(filename): print "file exists, exiting";sys.exit()
                uv = init_uv(filename)
            preamble = (p['uvw'],p['t'],p['bl'])
            uv['pol'] = a.miriad.str2pol[p['pol']]
            uv['lst'] = p['lst']
            uv['ra'] = p['lst']
            uv.write(preamble,d,f)
            cnt += 1    
        if n.sum(cont)<size: 
            del(uv)
            print "head node exiting..."
            break
#if rank==0:print len(datas),[len(d) for d in datas]
#if rank==0: print datas
#data = n.concatenate(comm.gather(data))
#if rank==0:
#    for p,d,f in data:
#        print p['t']
#
#
#data = n.sort(data,order='preamble')
#for p,d,f in data:
#    print p
#sys.exit()
#if rank==0:
#    cnt = 0
#    data = n.concatenate((datas,newdatas))
#    maxcnt=7e4
#    for p,d,f in data:
#        if cnt%maxcnt==0:
#            if cnt>0:del(uv)
#            filename = '.'.join(('sim',str(p['t']),'uv'))
#            if os.path.exists(filename): sys.exit()
#            uv = init_uv(filename)
#        preamble = (p['uvw'],p['t'],p['bl'])
#        uv['pol'] = a.miriad.str2pol[p['pol']]
#        uv['lst'] = p['lst']
#        uv['ra'] = p['lst']
#        uv.write(preamble,d,f)
#        cnt += 1
#    del(uv)
    
    
    
    