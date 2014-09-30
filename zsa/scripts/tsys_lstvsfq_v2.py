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
times, data, flags, stat = C.zsa.get_dict_of_uv_data(args, antstr, polstr, stats=stats, verbose=True)

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

bl_key = data.keys()[0] #to get size of array
pol_key = data[bl_key].keys()[0]
if opts.median:
    avg_var = []
    n_ints = []
else:
    avg_var = n.zeros_like(stat['var'][bl_key][pol_key])
    n_ints = n.zeros_like(stat['var'][bl_key][pol_key])
cnt = 0
for bl in stat['var'].keys():
    for p in stat['var'][bl].keys():
        if opts.median:
            avg_var.append(stat['var'][bl][p])
            n_ints.append(stat['var'][bl][p])
        else:
            avg_var += stat['var'][bl][p]
            n_ints += stat['cnt'][bl][p]
            cnt+=1

if opts.median:
    avg_var = n.median(avg_var, axis=0)
    n_ints = n.median(n_ints, axis=0)

else:
    avg_var = avg_var/cnt
    n_ints = n_ints/cnt


rescale = n.sqrt(inttime*sdf*1e9*2)#2 for pol.
#No factor of 2 in 2kt. 

freqs = freqs*1e3
extent = (freqs[0], freqs[-1], t2, t1)


TSYS = (1./1000)*n.sqrt(avg_var)*rescale*jy2T #in kelvin
TSYS_NOISE = TSYS/n.sqrt(n_ints)

TSYS164 = TSYS[:,130]
P.plot(lsts*12./n.pi,TSYS164)
P.title('Tsys 164 MHz')
P.xlabel('LST (Hours)')
P.ylabel('Tsys (K)')
P.show()
#P.plot(freqs, n.sum(n_ints, axis=0))
#P.show()

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
