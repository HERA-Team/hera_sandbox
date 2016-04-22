#! /usr/bin/env python
import aipy as a
import capo.arp as arp
import capo.frf_conv as fringe
import capo.zsa as zsa
import numpy as n, pylab as p
import sys, os, optparse
from scipy.special import erf

def skew(cenwid, bins):
        return n.exp(-(bins-cenwid[0])**2/(2*cenwid[1]**2))*(1+erf(cenwid[2]*(bins-cenwid[0])/(n.sqrt(2)*cenwid[1]))) 

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('--outpath', action='store',
    help='Output path to write to.')
o.add_option('--frpad',default=1.0,type=float,help='make the fringe rate convolution longer by this factor (default 1.0)')
o.add_option('--alietal', action='store_true',
        help='Uses normalization for alietal frf,(default=False)')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
nants = uv['nants']
inttime = uv['inttime'] * 4 #integration time hack
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
pol = a.miriad.pol2str[uv['pol']]


#manually calculate the inttime so frf time bins match data
(uvw,t1,(i,j)),d = uv.read()
(uvw,t2,(i,j)),d = uv.read()
while t1 == t2:
    (uvw,t2,(i,j)),d = uv.read()
inttime = (t2-t1)* (3600*24)
del(uv)

#inttime=50.

#set channel to make frf
mychan=101

##use calculated inttime to generate correct frf bins
frbins = fringe.gen_frbins(inttime)
#frbins = n.arange( -.5/inttime+5e-5/2, .5/inttime,5e-5)
#DEFAULT_FRBINS = n.arange(-.01+5e-5/2,.01,5e-5) # Hz

#Get only the antennas of interest
sep2ij, blconj, bl2sep = zsa.grid2ij(aa.ant_layout)
print "Using normalization for old FRF: ", opts.alietal

print "Looking for baselines matching ", opts.ant
ants = [ b[0] for b in a.scripting.parse_ants(opts.ant, nants) ]
seps = [ bl2sep[b] for b in ants ]
seps = n.unique(seps)
print 'These are the spearations that we are going to use ', seps
print 'This is the channel we are using to build the frf: ',mychan
print 'Current inttime use for gen_frbins: ',inttime
#Get the fir filters for the separation used.
firs = {}
for sep in seps:
    c = 0 
    while c != -1:
        ij = map(int, sep2ij[sep].split(',')[c].split('_'))
        bl = a.miriad.ij2bl(*ij)
        if blconj[bl]: c+=1
        else: break
    frp, bins = fringe.aa_to_fr_profile(aa, ij, mychan, bins=frbins,frpad=opts.frpad)
    #frp, bins = fringe.aa_to_fr_profile(aa, ij, mychan) ## for default fr_bins
    #timebins, firs[sep] = fringe.frp_to_firs(frp, bins, aa.get_afreqs(), fq0=aa.get_afreqs()[100],mdl=skew,startprms=(.001,.001,-50),frpad=opts.frpad)
    timebins, firs[sep] = fringe.frp_to_firs(frp, bins, aa.get_afreqs(), fq0=aa.get_afreqs()[mychan], startprms=(.001*opts.frpad,.0001))
    #timebins, firs[sep] = fringe.frp_to_firs(frp, bins, aa.get_afreqs(), fq0=aa.get_afreqs()[mychan],frpad=opts.frpad, alietal=opts.alietal )
    
baselines = ''.join(sep2ij[sep] for sep in seps)
times, data, flags = zsa.get_dict_of_uv_data(args, baselines, pol, verbose=True)
lsts = [ aa.sidereal_time() for k in map(aa.set_jultime(), times) ]

_d = {}
_w = {}
for bl in data.keys():
    if not _d.has_key(bl): _d[bl],_w[bl] = {}, {}
    #get filter which is baseline dependent.
    sep = bl2sep[bl]
    fir = firs[sep]
    if blconj[bl]: fir = n.conj(fir)
    print map(int, a.miriad.bl2ij(bl)), sep, blconj[bl]
    for pol in data[bl].keys():
        if not _d[bl].has_key(pol): _d[bl][pol], _w[bl][pol] = {}, {}
        _d[bl][pol] = n.zeros_like(data[bl][pol])
        _w[bl][pol] = n.zeros_like(data[bl][pol])
        for ch in xrange(data[bl][pol].shape[1]):
            flg = n.logical_not(flags[bl][pol][:,ch])
            #_d[bl][pol][:,ch] = n.convolve(flags[bl][pol][:,ch]*data[bl][pol][:,ch], n.conj(firs[ch,:]), mode='same')
            _d[bl][pol][:,ch] = n.convolve(flg*data[bl][pol][:,ch], n.conj(fir[ch,:]), mode='same')
            #_d[bl][pol][:,ch] = n.convolve(flg*data[bl][pol][:,ch], firs[ch,:], mode='same')
            #_w[bl][pol][:,ch] = n.convolve(flags[bl][pol][:,ch], n.abs(n.conj(firs[ch,:])), mode='same')
            _w[bl][pol][:,ch] = n.convolve(flg, n.abs(n.conj(fir[ch,:])), mode='same')
            #_w[bl][pol][:,ch] = n.convolve(flg, n.abs(firs[ch,:]), mode='same')
            #_d[bl][pol][:,ch] = n.where(flags[bl][pol][:,ch]>0, _d[bl][pol][:,ch]/_w[bl][pol][:,ch], 1)  
            _d[bl][pol][:,ch] = n.where(flg>0, _d[bl][pol][:,ch]/_w[bl][pol][:,ch], 1)  

def mfunc(uv, p, d, f):
    uvw,t,(i,j) = p
    index = n.where(times==t)
    bl = a.miriad.ij2bl(i,j)
    pol = a.miriad.pol2str[uv['pol']]
    try:
        #The flags are the same flags as input.
        d_ = _d[bl][pol][index,:].reshape(203)
    except(KeyError): return p,None,None
    else: return p, d_, f
        

for filename in args:
    outfile = filename+'L'
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, ants=baselines)
    if opts.outpath:
        uvo = a.miriad.UV(opts.outpath+'/'+outfile, status='new')
        print 'Writing %s'%(opts.outpath+'/'+outfile)
    else:
        uvo = a.miriad.UV(outfile, status='new')
        print 'Writing %s'%(outfile)
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+' inttime={0:.3f}s'.format(inttime)+' \\n', raw=True)
