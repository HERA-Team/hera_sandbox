#! /usr/bin/env python
import capo as C, pylab as P
import numpy as n
import aipy as a
import sys, optparse, capo

o = optparse.OptionParser()
o.add_option('--npz', action='store_true',
    help='read npz file instead of reading all the uv files again')
o.add_option('--median', action='store_true',   
    help='take median instead of average')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa('psa6240_v003', uv['sdf'], uv['sfreq'], uv['nchan'])
print uv['sdf']
freqs = aa.get_afreqs()
jy2T = capo.pspec.jy2T(freqs)
nchans = uv['nchan']

def get_dict_of_uv_data(filenames, antstr, polstr, decimate=1, decphs=0, verbose=False, recast_as_array=True, stat=None):
    times, dat, flg, stats = [], {}, {}, {}
    for st in stat:
        stats[st] = {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        if decimate > 1: uv.select('decimate', decimate, decphs)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            if len(times) == 0 or t != times[-1]: times.append(t)
            bl = a.miriad.ij2bl(i,j)
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            pol = a.miriad.pol2str[uv['pol']]
            if not dat[bl].has_key(pol):
                dat[bl][pol],flg[bl][pol] = [],[]
            dat[bl][pol].append(d)
            flg[bl][pol].append(f)
            for st in stat:
                if not stats[st].has_key(bl): stats[st][bl] = {}
                if not stats[st][bl].has_key(pol): stats[st][bl][pol] = []
                stats[st][bl][pol].append(uv[st])
    if recast_as_array:
        # This option helps reduce memory footprint, but it shouldn't
        # be necessary: the replace below should free RAM as quickly
        # as it is allocated.  Unfortunately, it doesn't seem to...
        for bl in dat.keys():
          for pol in dat[bl].keys():
            dat[bl][pol] = n.array(dat[bl][pol])
            flg[bl][pol] = n.array(flg[bl][pol])
            for st in stat:
                stats[st][bl][pol] = n.array(stats[st][bl][pol])
    return n.array(times), dat, flg, stats


sep='0,1'
antstr=''
antstr += capo.dfm.grid2ij(aa.ant_layout)[0][sep]
antstr = '1_4'
polstr = a.miriad.pol2str[uv['pol']]
inttime = uv['inttime']*4
sdf = uv['sdf']
del(uv)
nbls = len(antstr.split(','))

stats = ['var', 'cnt']

print 'antstr=', antstr
print 'polstr=', polstr
print 'inttime=', inttime
times, data, flags, stat = get_dict_of_uv_data(args, antstr, polstr, stat=stats, verbose=True)

import IPython 
IPython.embed()
lsts = []
for t in times:
    aa.set_jultime(t)
    lsts.append(aa.sidereal_time()) 
lsts = n.array(lsts)

#for plotting the correct lsts.
step = lsts[1] - lsts[0]
if lsts[0] > lsts[-1]:#wrapped
    diff = 2*n.pi - lsts[0]
    lstsmod = ((lsts + diff)%(2*n.pi)) - diff
    t1,t2 = (lstsmod[0]-0.5*step)*12/n.pi, (lstsmod[-1]+0.5*step)*12/n.pi
else:
    t1 = (lsts[0]-.5*step)*12/n.pi
    t2 = (lsts[-1]+.5*step)*12/n.pi
print t1,t2

extent = (freqs[0], freqs[-1], t2, t1)
counts = None
for bl in data.keys():
    print n.all(counts == stat['cnt'][bl][polstr])
    counts = stat['cnt'][bl][polstr]
    print n.max(counts)
    C.arp.waterfall(counts, mode='lin', extent=extent)
    P.colorbar(shrink=.5)
    P.title('%s_%s'%a.miriad.bl2ij(bl))
    P.xlabel('Freqs')
    P.ylabel('LST')
    P.show()

C.zsa.waterfall(TSYS, mode='lin', mx=1000, drng=1000, extent=extent)
P.title('Tsys in K')
P.xlabel('Frequency (MHz)')
P.ylabel('LST (Hours)')
P.figure(2)
P.title('Integrations per LST bin')
P.xlabel('Frequency (MHz)')
P.ylabel('LST (Hours)')
C.zsa.waterfall(n_ints,mode='lin', extent=extent)
P.figure(3)
P.title('Integrations per LST bin')
P.xlabel('Frequency (MHz)')
C.zsa.waterfall(TSYS_NOISE,mode='lin', extent=extent)
P.figure(3)
P.ylabel('LST (Hours)')
P.show()
