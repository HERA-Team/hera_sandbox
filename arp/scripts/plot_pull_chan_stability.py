#! /usr/bin/env python
import numpy as n, pylab as p, capo as C, aipy as a
import os

LST_RES = 43.
STARTJD = 2456242
NDAYS = 134

def lstbin(lst, lst_res=LST_RES):
    lst_res = lst_res / a.const.sidereal_day * (2*n.pi)
    return C.pspec.bin2uv(C.pspec.uv2bin(0,0,lst,lst_res=lst_res),lst_res=lst_res)[-1]

lstbins = n.arange(0,2*n.pi, LST_RES / a.const.sidereal_day * (2*n.pi))

#colors = 'kbmcrg'

densityfile = 'stability.npz'

if not os.path.exists(densityfile):
    import aipy as a, glob, sys
    ks = None
    ks = map(lambda t: str(a.miriad.ij2bl(*t)), [(1,4)])
    aa = a.cal.get_aa('psa6240_v003', n.array([.15]))
    sun = a.cal.get_catalog(srcs=['Sun'])['Sun']
    files = sys.argv[1:]
    files = glob.glob('zen.2456*.npz')
    density = n.zeros((lstbins.size, NDAYS), dtype=n.complex)
    wgt = n.zeros((lstbins.size, NDAYS), dtype=n.int)
    for filename in files:
        jd = float('.'.join(filename.split('.')[1:3]))
        aa.set_jultime(jd)
        sun.compute(aa)
        if sun.alt > -5 * a.ephem.degree:
            print 'Skipping', filename
            continue
        else:
            print 'Reading', filename
        day = int(jd) - STARTJD
        npz = n.load(filename)
        if ks is None: ks = npz.files
        for k in ks:
            if k.startswith('t'): continue
            try: lst,dat = npz['t'+k], npz[k]
            except(KeyError): continue
            #lst = n.where(lst > 5, lst-2*n.pi, lst)
            for L,D in zip(lst,dat):
                b = n.argmin(n.abs(lstbins - L))
                density[b,day] += D
                wgt[b,day] += 1
            #c = colors[int(k)%len(colors)]
            #x = n.around((lst+1)/(2*n.pi)*density.shape[0]).astype(n.int)
            #y = n.around(npz[k].real / 2 + density.shape[1]/2).astype(n.int)
            #density[x,y] += 1
            #print density.shape, density.max()
            #p.plot(lst, npz[k], c+',', alpha=.5)
            #p.plot(lst, npz[k], c+'-', linewidth=1, alpha=.5)
    #C.arp.waterfall(density, drng=1)
    #p.show()
    density = n.transpose(density)
    wgt = n.transpose(wgt)
    n.savez(densityfile, stability=density, wgt=wgt, lsts=lstbins, startjd=STARTJD)

npz = n.load(densityfile)
print npz.files
dat = npz['stability'][:,:2000]
wgt = npz['wgt']
#wgt = n.transpose(wgt) # XXX

_dat = n.fft.ifft(dat, axis=0)
_wgt = n.fft.ifft(wgt, axis=0)

for i in xrange(_dat.shape[1]):
    if i % 10 == 0: print i
    gain = n.sqrt(n.average(wgt[:,i]**2))
    if gain < .1: continue
    _dat[:,i],info = a.deconv.clean(_dat[:,i],_wgt[:,i],tol=1e-4)
    _dat[:,i] += info['res'] / gain

var = n.abs(_dat)**2
print var.shape
var.shape = (var.shape[0], var.shape[1]/25, 25)
var = n.average(var, axis=2)
print var.shape
_days = n.fft.fftfreq(var.shape[0])

p.figure(figsize=(5,4))
if False:
    for i in xrange(16,32):
        dat = var[:_days.size/2,i]
        #dat /= dat[0]
        p.loglog(_days[:_days.size/2]+.005, dat)
else:
    dat = n.average(var[:_days.size/2,16:32], axis=1)
    dat /= dat[0]
    p.loglog(_days[:_days.size/2]+.005, dat)
p.ylim(3e-5,3)
p.xlim(3e-3,1)
p.ylabel(r'Power [${\rm K}^2$, arbitrary scale]')
p.xlabel(r'Temporal Frequency [1/days]')
p.grid()
p.subplots_adjust(right=.95, top=.95, left=.17, bottom=.15)
p.show()

#p.figure(figsize=(12,6))
#p.subplot(211); C.arp.waterfall(dat, mode='log', drng=6)#, extent=rng, origin='lower', cmap='gist_earth_r', mode='lin', mx=40)
#p.subplot(212); C.arp.waterfall(var, mode='log', drng=6)#, extent=rng, origin='lower', cmap='gist_earth_r', mode='lin', mx=40)
##p.colorbar(ticks=[0,10,20,30,40])
##p.ylim(-500,500)
##p.ylim(-4500,4500)
##p.xlim(0,n.pi)
##p.xlim(0,2)
##p.grid()
##p.subplots_adjust(right=.99, bottom=.2)
##p.xlabel('LST [radians]')
##p.ylabel('Visibility (real) [Jy]')
#p.show()
