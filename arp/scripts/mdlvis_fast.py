#! /usr/bin/env python
'''Proof of concept for a fast visibility simulator that takes advantage of equal-spaced freq channels'''
import aipy as a, numpy as n, optparse, sys
import pylab as P
from aipy._cephes import i0

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('--nchan', dest='nchan', default=1024, type='int',
    help='Number of channels in simulated data if no input data to mimic.  Default is 1024')
o.add_option('--sfreq', dest='sfreq', default=.100, type='float',
    help='Start frequency (GHz) in simulated data if no input data to mimic.  Default is 0.100')
o.add_option('--sdf', dest='sdf', default=.100/1024, type='float',
    help='Channel spacing (GHz) in simulated data if no input data to mimic.  Default is .100/1024')
o.add_option('--inttime', dest='inttime', default=10, type='float',
    help='Integration time (s) in simulated data if no input data to mimic.  Default is 10')
o.add_option('--startjd', dest='startjd', default=2455600., type='float',
    help='Julian Date to start observation if no input data to mimic.  Default is 2454600')
o.add_option('--endjd', dest='endjd', default=2455601., type='float',
    help='Julian Date to end observation if no input data to mimic.  Default is 2454601')
opts,args = o.parse_args(sys.argv[1:])

NCHAN = opts.nchan
kaiser3 = lambda x: i0(n.pi * 3 * n.sqrt(1-(2*x/(NCHAN-1) - 1)**2)) / i0(n.pi * 3)
times = n.arange(opts.startjd, opts.endjd, opts.inttime/a.const.s_per_day)

aa = a.cal.get_aa(opts.cal, opts.sdf, opts.sfreq, opts.nchan)
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
pols = ['yy']
afreqs = aa.get_afreqs()
dlys = n.fft.fftfreq(NCHAN, d=afreqs[1]-afreqs[0])

def recenter(d, cen): return n.concatenate([d[cen:], d[:cen]])
def sim(tau): return n.exp(-2j*n.pi*afreqs*tau)

w = n.fromfunction(kaiser3, (opts.nchan,))
FOOTPRINT = 30 
SUBSAMPLE = 64

if True: # Calculate kernel explicitly
    _w = []
    #for tau in n.arange(-FOOTPRINT,FOOTPRINT+1,1./SUBSAMPLE):
    subx = []
    for tau in n.concatenate([n.arange(0,FOOTPRINT+1./SUBSAMPLE,1./SUBSAMPLE),n.arange(-FOOTPRINT,0,1./SUBSAMPLE)]):
        tau *= dlys[1]
        _w.append(n.fft.ifft(w*sim(tau))[0])
        subx.append(tau)
    _w = n.array(_w)
    subx = n.array(subx)
    #P.plot(n.abs(_w))
    #P.plot(n.angle(_w))
    #P.show()
else:
    mw = w
    if False:
        _w = n.fft.ifft(mw)
        subx = n.arange(0,NCHAN)
    else:
        _w = n.zeros(w.size * SUBSAMPLE, dtype=n.complex64)
        _w[:w.size] = mw
        #_w[:w.size/2] = w[w.size/2:]
        #_w[-w.size/2:] = w[:w.size/2]
        _w = n.fft.ifft(_w) * SUBSAMPLE
        subx = n.arange(0,NCHAN, 1./SUBSAMPLE)

taubin = 10
tau = dlys[taubin] + dlys[1]/8
# The default implementation
d = sim(tau)
_d = n.fft.ifft(d*w)
__d = n.fft.fft(_d) / n.where(w == 0, 1, w)
d__d = d - __d
# A gridded implementation
_g = n.zeros(NCHAN, dtype=n.complex64)
taubin = n.around(tau / dlys[1])
taures = tau / dlys[1] - taubin
#phs = n.exp(-2j*n.pi*taures) # Account for residual phase within bin
#_g[taubin] += phs
print int(taures*SUBSAMPLE)
for cnt in range(-FOOTPRINT,FOOTPRINT+1):
    #_g[taubin+cnt] += _w[cnt]
    _g[taubin+cnt] += _w[-cnt*SUBSAMPLE+int(n.around(taures*SUBSAMPLE))]
__g = n.fft.fft(_g) / n.where(w == 0, 1, w)
d__g = d - __g

P.subplot(411)
P.plot(n.abs(d), 'k-'); P.plot(d.real, 'k-.')#; P.plot(d.imag, 'k:')
P.xlim(0,NCHAN)
P.ylim(-1.1,1.1)

P.subplot(412)
_d /= n.abs(_d).max()
_g /= n.abs(_g).max()
_w /= n.abs(_w).max()#; _w = recenter(_w, -n.around(tau / dlys[1] * SUBSAMPLE))
P.plot(dlys, n.abs(_d), 'k-'); P.plot(dlys, _d.real, 'k-.'); P.plot(dlys, _d.imag, 'k:')
P.plot(dlys, n.abs(_g), 'r-'); P.plot(dlys, _g.real, 'r-.'); P.plot(dlys, _g.imag, 'r:')
P.plot(subx, n.abs(_w), 'g-'); P.plot(subx, _w.real, 'g-.'); P.plot(subx, _w.imag, 'g:')
P.xlim(0,2*taubin*dlys[1])

P.subplot(413)
P.plot(n.abs(__d), 'k-'); P.plot(__d.real, 'k-.')#; P.plot(__d.imag, 'k:')
P.plot(n.abs(__g), 'r-'); P.plot(__g.real, 'r-.')#; P.plot(__g.imag, 'r:')
P.xlim(0,NCHAN)
P.ylim(-1.1,1.1)

P.subplot(414)
P.semilogy(n.abs(d__d), 'k-')
P.semilogy(n.abs(d__g), 'r-')
P.xlim(0,NCHAN)
P.ylim(1e-5, 2e0)
P.grid()

P.show()
    

F = 3
curtime = None
L = 33
#import capo as C; C.pfb.__set_pm__(NCHAN*F, 'kaiser3',1,1.); window = C.pfb.pm['window']
#window = n.fromfunction(kaiser3, (NCHAN,))
w = n.fromfunction(kaiser3, (L,))
#w = n.ones_like(w)
print w.sum()
#w /= w.sum() * 8 / 1.1 * 1.02
SUBSAMPLE = F * 10
FOOTPRINT=10
_w = n.zeros(w.size * SUBSAMPLE, dtype=n.complex64)
_w[:w.size/2+1] = w[w.size/2:]
_w[-w.size/2+1:] = w[:w.size/2]
_w = n.fft.ifft(_w)
norm = n.zeros(NCHAN*F, dtype=n.complex64)
for cnt in range(-FOOTPRINT,FOOTPRINT+1):
    norm[cnt] += _w[SUBSAMPLE*cnt]
norm = n.fft.fft(norm)
if False:
    P.subplot(211)
    P.semilogy(n.abs(_w))
    #P.plot(_w.real)
    #P.plot(_w.imag)
    P.subplot(212)
    thing = n.zeros(w.size * SUBSAMPLE, dtype=n.complex64); thing[100] = 1.
    P.plot(n.abs(thing))
    P.plot(n.abs(n.fft.ifft(n.fft.fft(thing) * _w)))
    P.show()
dlys = n.fft.fftfreq(int(NCHAN*F), d=afreqs[1]-afreqs[0])
def mdl(uv, p, d, f):
    global curtime, eqs
    uvw, t, (i,j) = p
    if i == j: return p, d, f
    if curtime != t:
        curtime = t
        aa.set_jultime(t)
        cat.compute(aa)
        eqs = cat.get_crds('eq', ncrd=3)
        flx = n.ones_like(cat.get_jys())
        aa.sim_cache(eqs, flx)
    #x,y,tau_g = aa.get_baseline(i,j, src=cat.values()[0])
    #tau_g = 122./3
    tau_g = 1
    tau_e,off = n.array(aa[j]._phsoff) - n.array(aa[i]._phsoff)
    taubin = (-tau_g + tau_e) / dlys[1]
    tauind = n.around(taubin)
    taubin -= tauind
    _d = n.zeros(NCHAN*F, dtype=n.complex64)
    phs = 1
    #phs = -n.exp(-1j*n.pi*taubin) # Account for residual phase within bin
    for cnt in range(-FOOTPRINT,FOOTPRINT+1):
        _d[tauind+cnt] += phs * _w[int(SUBSAMPLE*(cnt+taubin))]
    #_d[tauind+1] += phs
    #_d[tauind-1] += phs
    d = n.fft.fft(_d)
    d /= norm
    #d /= d.max()
    d = n.concatenate([d[-d.size/2:],d[:d.size/2]]) 
    #d = d[int(NCHAN*(F-1)/2):int(NCHAN*(F+1)/2)]# / partwin
    print d.shape
    sd = n.zeros(int(NCHAN*F), dtype=n.complex64)
    #sd[int(NCHAN*(F-1)/2):int(NCHAN*(F+1)/2)] = aa.sim(i, j, pol=a.miriad.pol2str[uv['pol']])
    sd[int(NCHAN*(F-1)/2):int(NCHAN*(F+1)/2)] = n.exp(-2j*n.pi*afreqs*(-tau_g)) #* window
    #sd = n.exp(-2j*n.pi*afreqs*(-tau_g))
    #asd = n.abs(sd)
    #sd /= n.where(asd == 0, 1, asd)
    _sd = n.fft.ifft(sd)
    #sd = sd[int(NCHAN*(F-1)/2):int(NCHAN*(F+1)/2)]
    print sd.shape
    print tau_g + tau_e, dlys[1], (tauind, taubin), _d[tauind], _sd[tauind] / n.abs(_sd[tauind])
    #print n.angle(sd * n.conj(d))
    #d = n.exp(2j*n.pi*afreqs*(tau_g + tau_e))
    P.subplot(311)
    P.plot( n.angle(d), 'k')
    #P.plot(afreqs, d.real, 'k')
    #P.plot(afreqs, d.imag, 'k:')
    P.plot( n.angle(sd), 'r')
    #P.plot( n.angle(sd*n.conj(d)), 'g')
    #P.plot(afreqs, sd.real, 'r')
    #P.plot(afreqs, sd.imag, 'r:')
    P.subplot(312)
    P.plot(n.abs(d), 'k')
    #P.plot(afreqs, d.real, 'k')
    #P.plot(afreqs, d.imag, 'k:')
    P.plot(n.abs(sd), 'r')
    #P.plot(afreqs, sd.real, 'r')
    #P.plot(afreqs, sd.imag, 'r:')
    P.subplot(313)
    #P.plot(dlys, n.abs(_d), 'k')
    #P.plot(dlys, n.abs(_sd), 'r')
    P.semilogy(dlys, n.abs(_d), 'k.-')
    P.semilogy(dlys, n.abs(_sd), 'r.-')
    #P.semilogy(dlys0, n.abs(n.fft.ifft(d)), 'g.-')
    P.xlim(-100,100)
    P.show()
    if not opts.flag: f = no_flags
    return p, n.where(f, 0, d), f



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
antpos = n.array([ant.pos for ant in aa], dtype=n.double)
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

