#! /usr/bin/env python
'''Proof of concept for a fast visibility simulator that takes advantage of equal-spaced freq channels'''
import aipy, numpy as np, optparse, sys
#import pylab as P

o = optparse.OptionParser()
aipy.scripting.add_standard_options(o, cal=True)
o.add_option('--freq', dest='freq', default=.150, type='float',
    help='Frequency (GHz) of map.  Default is 0.15')
o.add_option('--inttime', dest='inttime', default=10, type='float',
    help='Integration time (s) in simulated data if no input data to mimic.  Default is 10')
o.add_option('--startjd', dest='startjd', default=2455600., type='float',
    help='Julian Date to start observation if no input data to mimic.  Default is 2454600')
o.add_option('--endjd', dest='endjd', default=2455601., type='float',
    help='Julian Date to end observation if no input data to mimic.  Default is 2454601')
opts,args = o.parse_args(sys.argv[1:])

NANT = 32
FREQ = opts.freq
START_JD = opts.startjd
END_JD = opts.endjd
INT_TIME = opts.inttime
POL = 'yy'

times = np.arange(START_JD, END_JD, INT_TIME / aipy.const.s_per_day)
aa = aipy.cal.get_aa(opts.cal, np.array([FREQ]))
pols = [POL]
aa.set_active_pol(pols[0])
afreqs = aa.get_afreqs()
antpos = {i:aa.get_baseline(0,i,'z').astype(np.float32) for i in xrange(NANT)}
bls = [(i,j) for i in antpos for j in antpos if i <= j]

hpm = aipy.map.Map(fromfits=sys.argv[-1]) # Map in eq coord, one freq?
px = np.arange(hpm.npix())
I_s = np.sqrt(hpm[px]) # sqrt so that (AI)^2 has one factor of sky
crd_equ = np.array(hpm.px2crd(px))

data = {(i,j,pol): [] for i,j in bls for pol in pols}

import time
print '# Antennas:', len(antpos)
print '# Baselines:', len(data)
print 'NSIDE:', hpm.nside()
print 'Starting', time.time()
for ti,jd in enumerate(times):
    print ti,'/',len(times)
    aa.set_jultime(jd)
    tx,ty,tz = crd_top = np.dot(aa.eq2top_m, crd_equ).astype(np.float32)
    crd_top = crd_top[:,tz > 0]
    # XXX assuming all beams are same for now
    A_s = aa[0].bm_response(crd_top, pol=pols[0][0]) # XXX fix pol treatment
    # XXX deal with antenna gains/phs
    AI_s = A_s * I_s[tz > 0]
    AI_s.shape = (1,-1)
    factor = np.array(-2j*np.pi*FREQ / aipy.const.c,dtype=np.complex64)
    antvis = np.empty((len(antpos),AI_s.size), dtype=np.complex64)
    for i in antpos: np.dot(factor*antpos[i], crd_top, out=antvis[i])
    np.exp(antvis, out=antvis)
    antvis *= AI_s
    data = np.zeros((len(antpos),len(antpos)), dtype=np.complex64)
    for i in antpos: np.dot(antvis[i:i+1].conj(), antvis[i:].T, out=data[i:i+1,i:])
    np.conj(data, out=data)
        
print 'Done', time.time()
sys.exit(0)

uv = a.miriad.UV('new.uv', status='new')
uv._wrhd('obstype','mixed-auto-cross')
uv._wrhd('history','MDLVIS: ' + ' '.join(sys.argv) + '\n') 
uv.add_var('telescop','a'); uv['telescop'] = 'AIPY'
uv.add_var('operator','a'); uv['operator'] = 'AIPY'
uv.add_var('version' ,'a'); uv['version'] = '0.0.2'
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
antpos = np.array([ant.pos for ant in aa], dtype=n.double)
uv['antpos'] = antpos.transpose().flatten()
uv.add_var('sfreq'   ,'d'); uv['sfreq'] = opts.sfreq
uv.add_var('freq'    ,'d'); uv['freq'] = opts.sfreq
uv.add_var('restfreq','d'); uv['restfreq'] = opts.sfreq
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

# Now start generating data
no_data = n.zeros(opts.nchan, dtype=n.complex64)
no_flags = n.zeros(opts.nchan, dtype=n.int32)
times = n.arange(opts.startjd, opts.endjd, opts.inttime/a.const.s_per_day)
for cnt,t in enumerate(times):
    print 'Timestep %d / %d' % (cnt+1, len(times))
    aa.set_jultime(t)
    uv['lst'] = aa.sidereal_time()
    uv['ra'] = aa.sidereal_time()
    uv['obsra'] = aa.sidereal_time()
    for i,ai in enumerate(aa):
        for j,aj in enumerate(aa):
            if j < i: continue
            crd = ai.pos - aj.pos
            preamble = (crd, t, (i,j))
            for pol in pols:
                uv['pol'] = a.miriad.str2pol[pol]
                preamble,data,flags = mdl(uv, preamble, None, None)
                if data is None:
                    data = no_data
                    flags = no_flags
                uv.write(preamble, data, flags)
del(uv)

