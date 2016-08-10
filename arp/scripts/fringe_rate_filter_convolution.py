#! /usr/bin/env python
'''Fringe Rate Filter Hack'''
import aipy as a, numpy as n, sys, os, optparse
import capo as C
import pylab as p

def grid2ij(GRID):
    '''
        bl_str = given sep, returns bls in string format.
        bl_conj = given a baseline (miriad bl) gives separation.
        bl2sep_str = given baseline (miriad) return its separation.    
    '''
    bls, conj = {}, {}
    for ri in range(GRID.shape[0]):
        for ci in range(GRID.shape[1]):
            for rj in range(GRID.shape[0]):
                for cj in range(GRID.shape[1]):
                    if ci > cj: continue
#                    if ri > rj and ci == cj: continue
#                    if ci > cj and ri == rj: continue
                    sep = (rj-ri, cj-ci)
                    sep = '%d,%d'%sep
                    i,j = GRID[ri, ci], GRID[rj,cj]
                    bls[sep] = bls.get(sep,[]) + [(i,j)]
    for sep in bls.keys():
        if sep == '0,0' or len(bls[sep]) < 2 or (sep[-1] == '0' and sep[0] == '-'): del(bls[sep]) 
    for sep in bls:
        conj[sep] = [i>j for i,j in bls[sep]]

    bl_str,bl_conj,bl2sep_str = {}, {}, {}
    for sep in bls:
        bl_str[sep],bl_list = [], []
        for (i,j),c in zip(bls[sep],conj[sep]):
            if c: i,j = j,i
            bl_list.append(a.miriad.ij2bl(i,j))
            bl_str[sep].append('%d_%d'%(i,j))
            bl2sep_str[a.miriad.ij2bl(i,j)] = bl2sep_str.get(a.miriad.ij2bl(i,j),'') + sep
            bl_conj[a.miriad.ij2bl(i,j)] = c
        bls[sep] = bl_list
        bl_str[sep] = ','.join(bl_str[sep])
    return bl_str,bl_conj,bl2sep_str


o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)#, ant=True, pol=True)
o.add_option('--sep', type='str', default='0,1',
    help='Separation type. Can give list delimitted by ;\n (rowdiff,coldiff) where col diff is positive\
          and (0,0) is upper left of array.')
o.add_option('-w', '--filt_width', type='int', default=401,
    help='Filter width in time domain in number of integration units')
o.add_option('--outpath', action='store', default='',
    help='output path.')

opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
print "inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly",
inttime = uv['inttime'] * 4 # XXX this is a hack for *E files that have inttime set incorrectly
print inttime
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
pol = a.miriad.pol2str[uv['pol']]
del(uv)

#Get only the antennas of interest
bl_str, conj, bl2sep = grid2ij(aa.ant_layout)
separations = opts.sep.split(';')
#get appropriate baseline to compute filter.
seps = {}
for sep in separations:
    bs = bl_str[sep].split(',')
    for sbl in bs:
        i,j = map(int, sbl.split('_'))
        if conj[a.miriad.ij2bl(i,j)]:continue
        else:
            seps[sep] = (i,j)
            break
        
baselines = ''
for sep in separations:
    baselines += bl_str[sep]
print baselines
print seps
    
#seps = {'0,1':(1,4), '1,1':(1,18), '-1,1':(1,48)}

#Get filters with proper conjugation etc. Using one for each baseline type
#too much time to generate al lof them.
filters = {}
frspace_filters = {}
for sep in seps.keys():
    beam_w_fr = C.frf_conv.get_beam_w_fr(aa, seps[sep])
    t, firs, frbins,frspace = C.frf_conv.get_fringe_rate_kernels(beam_w_fr, inttime, opts.filt_width)
    filters[sep] = firs
    frspace_filters[sep] = frspace

print filters

#
_d = {}
_w = {}
times, data, flags = C.arp.get_dict_of_uv_data(args, baselines, pol, verbose=True)
lsts = []
for t in times:
    aa.set_jultime(t)
    lsts.append(aa.sidereal_time())

#print lsts
#ii,jj = seps[seps.keys()[0]]
#testbl = a.miriad.ij2bl(ii, jj)
#print data[testbl]['I'].shape
#C.arp.waterfall(data[testbl]['I'], extent=(0,203,lsts[-1],lsts[0]));p.colorbar(shrink=.5)
#p.show()
#frates = n.fft.fftshift(n.fft.fftfreq(len(data[testbl]['I']), 42.8) * 1e3)
#C.arp.waterfall(n.fft.fftshift(n.fft.ifft(n.fft.fftshift(data[testbl]['I'],axes=0),axis=0),axes=0),extent=(0,203,frates[-1],frates[1]));p.colorbar(shrink=.5)
#p.show()


for bl in data.keys():
    if not _d.has_key(bl): _d[bl],_w[bl] = {}, {}
    #get filter which is baseline dependent.
    b1,b2 = map(int,a.miriad.bl2ij(bl))
    sep = bl2sep[bl]
    firs = filters[sep]
    if conj[bl]: firs = n.conj(firs)
    firs = n.conj(firs) # XXX check where this sign change came from
    print b1,b2, sep, conj[bl]
    for pol in data[bl].keys():
        if not _d[bl].has_key(pol): _d[bl][pol], _w[bl][pol] = {}, {}
        _d[bl][pol] = n.zeros_like(data[bl][pol])
        _w[bl][pol] = n.zeros_like(data[bl][pol])
        for ch in xrange(data[bl][pol].shape[1]):
            flg = n.logical_not(flags[bl][pol][:,ch])
            #_d[bl][pol][:,ch] = n.convolve(flags[bl][pol][:,ch]*data[bl][pol][:,ch], n.conj(firs[ch,:]), mode='same')
            _d[bl][pol][:,ch] = n.convolve(flg*data[bl][pol][:,ch], n.conj(firs[ch,:]), mode='same')
            #_d[bl][pol][:,ch] = n.convolve(flg*data[bl][pol][:,ch], firs[ch,:], mode='same')
            #_w[bl][pol][:,ch] = n.convolve(flags[bl][pol][:,ch], n.abs(n.conj(firs[ch,:])), mode='same')
            _w[bl][pol][:,ch] = n.convolve(flg, n.abs(n.conj(firs[ch,:])), mode='same')
            #_w[bl][pol][:,ch] = n.convolve(flg, n.abs(firs[ch,:]), mode='same')
            #_d[bl][pol][:,ch] = n.where(flags[bl][pol][:,ch]>0, _d[bl][pol][:,ch]/_w[bl][pol][:,ch], 1)  
            # XXX how do deal with uneven time-domain sampling?
            #_d[bl][pol][:,ch] = n.where(flg>0, _d[bl][pol][:,ch]/_w[bl][pol][:,ch], 0)  

#C.arp.waterfall(n.fft.fftshift(n.fft.ifft(n.fft.fftshift(_d[testbl]['I'],axes=0),axis=0),axes=0),extent=(0,203,frates[-1],frates[1]));p.colorbar(shrink=.5)
#p.show()

def mfunc(uv, p, d, f):
    uvw,t,(i,j) = p
    index = n.where(times==t)
    bl = a.miriad.ij2bl(i,j)
    pol = a.miriad.pol2str[uv['pol']]
    try:
        #The flags are the same flags as input.
        d_ = _d[bl][pol][index,:].reshape(d.size)
    except(KeyError): return p,None,None
    else: return p, d_, f
        

for filename in args:
    outfile = filename+'L'
    print 'Writing %s'%(opts.outpath+'/'+outfile)
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, ants=baselines)
    if opts.outpath:
        uvo = a.miriad.UV(opts.outpath+'/'+outfile, status='new')
    else:
        uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+'\n', raw=True)





