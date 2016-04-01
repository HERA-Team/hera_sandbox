#! /usr/bin/env python
import capo.hex as hx, capo.arp as arp, capo.red as red, capo.omni as omni
import numpy as n, pylab as p, aipy as a
import sys,optparse

o = optparse.OptionParser()
o.add_option('--cal', action='store',
    help='File path for the connections.')
o.add_option('--plot', action='store_true', help='Plot things.')
opts,args = o.parse_args(sys.argv[1:])
connection_file=opts.cal
PLOT=opts.plot

def flatten_reds(reds):
    freds = []
    for r in reds:
        freds += r
    return freds

def save_sols(s):
    s2 = {}
    for k,i in s.iteritems():
        s2[str(k)] = i
    n.savez('fcsols.npz',**s2)

#hera info assuming a hex of 19 and 128 antennas
info = hx.hera_to_info(3, 128, connections=connection_file)
infotest = hx.hera_to_info(3, 128, connections=connection_file, ubls=[(80,104)])
reds = flatten_reds(info.get_reds())
redstest = infotest.get_reds()#for plotting 

#Read in data here.
ant_string =','.join(map(str,info.subsetant))
bl_string = ','.join(['_'.join(map(str,k)) for k in reds])
times, data, flags = arp.get_dict_of_uv_data(args, bl_string, 'xx', verbose=True)
dataxx = {}
for (i,j) in data.keys():
    dataxx[(i,j)] = data[(i,j)]['xx']
fqs = n.linspace(.1,.2,1024)

#gets phase solutions per frequency.
fc = omni.FirstCal(dataxx,fqs,info)
sols = fc.run()
save_sols(sols)
#save solutions
dataxx_c = {}
for (a1,a2) in info.bl_order():
    if (a1,a2) in dataxx.keys():
        dataxx_c[(a1,a2)] = dataxx[(a1,a2)]*omni.get_phase(fqs,sols[a1])*n.conj(omni.get_phase(fqs,sols[a2]))
    else:
        dataxx_c[(a1,a2)] = dataxx[(a2,a1)]*omni.get_phase(fqs,sols[a2])*n.conj(omni.get_phase(fqs,sols[a1]))

#def waterfall(d, ax, mode='log', mx=None, drng=None, recenter=False, **kwargs):
#    if n.ma.isMaskedArray(d): d = d.filled(0)
#    if recenter: d = a.img.recenter(d, n.array(d.shape)/2)
#    d = arp.data_mode(d, mode=mode)
#    if mx is None: mx = d.max()
#    if drng is None: drng = mx - d.min()
#    mn = mx - drng
#    return ax.imshow(d, vmax=mx, vmin=mn, aspect='auto', interpolation='nearest', **kwargs)
#
#plotting data
redbls = []
for r in redstest: redbls += r
redbls = n.array(redbls)
#print redbls.shape
#dm = divmod(len(redbls), n.round(n.sqrt(len(redbls))))
#nr,nc = int(dm[0]),int(dm[0]+n.ceil(float(dm[1])/dm[0]))
#fig,ax = p.subplots(nrows=nr,ncols=nc,figsize=(14,10))
#for i,bl in enumerate(redbls):
#    bl = (bl[0],bl[1])
#    try: 
#        waterfall(dataxx[bl], ax[divmod(i,nc)], mode='phs')
#        ax[divmod(i,nc)].set_title('%d,%d'%(bl))
#    except(KeyError):
#        waterfall(dataxx[bl[::-1]], ax[divmod(i,nc)], mode='phs')
#        ax[divmod(i,nc)].set_title('%d,%d'%(bl[::-1]), color='m')
#fig.subplots_adjust(hspace=.5)

if PLOT:
    for bl in redbls:
        bl = tuple(bl)
        try:
            #p.subplot(211); arp.waterfall(dataxx[bl], mode='log',mx=0,drng=3); p.colorbar(shrink=.5)
            #p.subplot(212); arp.waterfall(dataxx_c[bl], mode='log',mx=0,drng=3); p.colorbar(shrink=.5)
            p.subplot(211); arp.waterfall(dataxx[bl], mode='phs'); p.colorbar(shrink=.5)
            p.subplot(212); arp.waterfall(dataxx_c[bl], mode='phs'); p.colorbar(shrink=.5)
            print sols[bl[0]] - sols[bl[1]]
            print bl
        except(KeyError):
            p.subplot(211); arp.waterfall(n.conj(dataxx[bl[::-1]]), mode='phs'); p.colorbar(shrink=.5)
            p.subplot(212); arp.waterfall(n.conj(dataxx_c[bl]), mode='phs'); p.colorbar(shrink=.5)
            print bl

        p.show()
import IPython; IPython.embed()
   





