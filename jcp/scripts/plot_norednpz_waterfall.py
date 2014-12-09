#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
import sys, optparse, os

class DataQuery:
    def __init__(self,fig,bls):
        self.x = 0
        self.y = 0
        self.cid = fig.canvas.mpl_connect('button_press_event',self)
    def __call__(self,event):
       try:
           self.y = event.ydata- 1
           line = n.floor(self.y)
           bl = bls[line]
           print a.miriad.bl2ij(bl)
       except(TypeError): pass

o = optparse.OptionParser()
o.add_option('--umax', type='float', default=n.Inf,
    help='Only show baselines shorter than this value')   
o.add_option('--umin', type='float', default=0.,
    help='Only show baselines longer than this value')   
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

fq = .16
aa = a.cal.get_aa(opts.cal, n.array([fq]))

magdict, wgtdict = {},{}
mags,bls,conj = [],[],[]
for npzfile in args:
    print 'Reading...', npzfile
    dat = n.load(npzfile)
    kpl = dat['kpl']
    keys = dat.files[:]
    keys.remove('kpl')

    for bl in keys:
        if bl[0] == 'w': continue
        wbl = 'w'+str(bl)
        i,j = a.miriad.bl2ij(bl)
        if i == 40 or j == 40: continue
        if i == 55 or j == 55: continue
        crd = aa.get_baseline(i,j)*fq
        mag = n.sqrt(crd[0]**2 + crd[1]**2)
        magdict[mag] = magdict.get(mag,0) + dat[bl]
        wgtdict[mag] = wgtdict.get(mag,0) + dat[wbl]
        mags.append(mag)
        bls.append(bl)
        conj.append(0)

mags,bls = n.unique(n.array(mags)),n.unique(n.array(bls))
inds = n.argsort(mags)
mags,bls,conj = n.take(mags,inds), n.take(bls,inds), n.take(conj,inds)
valid = n.where(mags > opts.umin, 1, 0) * n.where(mags < opts.umax, 1, 0)
mags,bls,conj = mags.compress(valid), bls.compress(valid), conj.compress(valid)

blmags = []
for bl in bls:
    i,j = a.miriad.bl2ij(bl)
    crd = aa.get_baseline(i,j)*fq
    mag = n.sqrt(crd[0]**2 + crd[1]**2)
    blmags.append(mag)

wfall = n.zeros((len(mags),len(kpl)))
for line,mag in enumerate(mags):
    spec = magdict[mag] / wgtdict[mag]
    inds = n.argsort(kpl)
    spec = n.take(spec,inds)
    wfall[line,:] = spec
    #if line == 0: p.semilogy(kpl[inds],spec);p.show()

kpl.sort()

fig = p.figure()
axis = fig.add_subplot(111)

p.imshow(n.log10(n.abs(wfall)),aspect='auto',interpolation='nearest',extent=(kpl[0],kpl[-1],len(mags),0),vmin=11,vmax=16)
p.colorbar()

query = DataQuery(fig,bls)

p.show()
