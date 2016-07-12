'''Fringe Rate Filter Hack'''
import frf_conv 
import aipy as a, numpy as n, sys, os, optparse
import capo as C,pylab as p

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('-w', '--filt_width', type='int', default=401,
    help='Filter width in time domain in number of integration units')

opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
print "inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly",
inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly
print inttime
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
pol = a.miriad.pol2str[uv['pol']]
del(uv)

#Get only the antennas of interest
sep2ij = {}
ij2sep = {}
s = ''
seps = {'0,1':(1,4), '1,1':(1,18), '-1,1':(1,48)}
for sep in seps.keys():
    sep2ij[sep] = C.dfm.grid2ij(aa.ant_layout)[0][sep].split(',')
    s += C.dfm.grid2ij(aa.ant_layout)[0][sep] + ','

    for ij in sep2ij[sep]:
        ij2sep[ij] = sep

conj = {}
toconj = C.dfm.grid2ij(aa.ant_layout)[1]
for k in toconj.keys():
    conj['%d_%d'%a.miriad.bl2ij(k)] = toconj[k]

#Get filters with proper conjugation etc. Using one for each baseline type
#too much time to generate al lof them.
filters = {}
for sep in seps.keys():
    beam_w_fr = frf_conv.get_beam_w_fr(aa, seps[sep])
    t, firs, frbins,frspace = frf_conv.get_fringe_rate_kernels(beam_w_fr, inttime, opts.filt_width)
    filters[sep] = firs

#
_d = {}
_w = {}
times, data, flags = C.arp.get_dict_of_uv_data(args, s, pol, verbose=True)

for bl in data.keys():
    if not _d.has_key(bl): _d[bl],_w[bl] = {}, {}
    #get filter which is baseline dependent.
    b1,b2 = map(int,a.miriad.bl2ij(bl))
    sep = ij2sep['%d_%d'%(b1,b2)]
    firs = filters[sep]
    print b1,b2, sep
    for pol in data[bl].keys():
        if not _d[bl].has_key(pol): _d[bl][pol], _w[bl][pol] = {}, {}
        _d[bl][pol] = n.zeros_like(data[bl][pol])
        _w[bl][pol] = n.zeros_like(data[bl][pol])
        for ch in xrange(data[bl][pol].shape[1]):
            flg = n.logical_not(flags[bl][pol][:,ch])
            #_d[bl][pol][:,ch] = n.convolve(flags[bl][pol][:,ch]*data[bl][pol][:,ch], n.conj(firs[ch,:]), mode='same')
            _d[bl][pol][:,ch] = n.convolve(flg*data[bl][pol][:,ch], n.conj(firs[ch,:]), mode='same')
            #_w[bl][pol][:,ch] = n.convolve(flags[bl][pol][:,ch], n.abs(n.conj(firs[ch,:])), mode='same')
            _w[bl][pol][:,ch] = n.convolve(flg, n.abs(n.conj(firs[ch,:])), mode='same')
            #_d[bl][pol][:,ch] = n.where(flags[bl][pol][:,ch]>0, _d[bl][pol][:,ch]/_w[bl][pol][:,ch], 1)  
            _d[bl][pol][:,ch] = n.where(flg>0, _d[bl][pol][:,ch]/_w[bl][pol][:,ch], 1)  
p.figure(1)
C.arp.waterfall(_d[517]['I'])
p.figure(2)
C.arp.waterfall(n.fft.ifft(_d[517]['I'],axis=0))
p.show()

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
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+'\n', raw=True)





