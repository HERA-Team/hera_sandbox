#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
import sys, optparse, os
from matplotlib import rc
rc('text',usetex=True)
rc('font', size=16)

def dk_du(z):
    '''2pi * [h Mpc^-1] / [wavelengths]'''
    return 2*n.pi / C.pspec.dL_dth(z)

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

fq = .1525
aa = a.cal.get_aa(opts.cal, n.array([fq]))
#print n.sqrt(aa.get_baseline(7,37)[0]**2 + aa.get_baseline(7,37)[1]**2)

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
        #if i == 7 and j == 37: print i,j, mag
        magdict[mag] = magdict.get(mag,0) + dat[bl]
        wgtdict[mag] = wgtdict.get(mag,0) + dat[wbl]
        mags.append(mag)
        bls.append(bl)
        if crd[0] < 0.:
            conj.append(1)
        else:
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
    #if i == 7 and j == 37: print i,j,mag
    blmags.append(mag)

#hi, bins = .14, 750 #published version
hi, bins = .14, 300
#hi, bins = .083, 600
step = hi/bins
kprs = n.arange(0,hi,step)
half = len(kpl)/2

wfall = n.zeros((len(kprs),half))
wgtfall = n.zeros((len(kprs),half))
hors,phors = n.zeros_like(kprs),n.zeros_like(kprs)
for ind,mag in enumerate(mags):
    spec = magdict[mag] / wgtdict[mag]
    #conjugate [doesnt work]
    if conj[ind] == 1: 
         spec = spec[::-1]
    inds = n.argsort(kpl)
    spec = n.take(spec,inds)
    kpr = dk_du(C.pspec.f2z(fq))*mag
    kpr = n.round(kpr/(step))*(step)
    #calculate the horizon
    hor = C.pspec.dk_deta(C.pspec.f2z(fq))*(mag/fq)
    phor = C.pspec.dk_deta(C.pspec.f2z(fq))*((mag/fq)+50)
    line = n.where(kprs == kpr)[0]
    hors[line] = hor
    phors[line] = phor
    #fold over kpl = 0
    spec = (spec[0:half]+spec[half:][::-1])/2 
    wfall[line,:] += spec
    wgtfall[line,:] += n.ones_like(spec)
    #if line == 0: p.semilogy(kpl[inds],spec);p.show()

wfall[wgtfall>0] /= wgtfall[wgtfall>0]
kpl.sort()

fig = p.figure()
axis = fig.add_subplot(111)

wfall = wfall.T
#multiply by 2 for beam cludge
#p.imshow(n.log10(2*n.abs(wfall)),aspect='auto',interpolation='nearest',extent=(0,kprs[-1],0,n.max(kpl)),vmin=11,vmax=20)
#p.imshow(n.log10(2*n.abs(wfall)),interpolation='nearest',extent=(0,kprs[-1],0,n.max(kpl)),vmin=11,vmax=15,cmap='gray')
p.imshow(n.log10(2*n.abs(wfall)),interpolation='nearest',extent=(0,kprs[-1],0,n.max(kpl)),vmin=11,vmax=15)
#cb = p.colorbar(orientation='horizontal',pad=.125,ticks=[11,13,15])
cb = p.colorbar(ticks=[11,13,15])
cb.set_label(r'${\rm log_{10}}[{\rm mK}^2 (h^{-1}{\rm Mpc})^3]$', size=12)

#save wedge data to npz
#n.savez('wedgedata.npz',data=n.log10(2*n.abs(wfall)),kprs=kprs,kpls=kpl)

#p.plot(kprs,hors,'.',lw=2,color='white')
p.scatter(kprs,hors,s=1,color='black')
#p.plot(kprs,phors,'.',lw=2,color='orange')

kcircs = n.arange(0.,0.5,0.05)
for kcirc in kcircs:
    kperps = n.linspace(0.,kcirc,1000.)
    kpls = []
    for kperp in kperps:
        if not kperp > kcirc:
            kpl = n.sqrt(kcirc**2 - kperp**2)
            kpls.append(kpl)
    kperps,kpls = n.array(kperps),n.array(kpls)
    p.plot(kperps,kpls,color='white')
    

query = DataQuery(fig,bls)

#p.title(r'$P(k)\ [{\rm mK}^2 (h^{-1}{\rm Mpc})^{3}]$')
p.title(r'${\rm log}_{10}[P(k)]$')
p.ylabel(r'$k_{\parallel}\ [h{\rm Mpc}^{-1}]$')
p.xlabel(r'$k_{\perp}\ [h{\rm Mpc}^{-1}]$')
#p.xticks([0,0.04,0.08,0.12])
p.xticks([0,0.06,0.12])
#p.ylim(0,1)
p.ylim(0,.5)
p.xlim(0,.12)

#ay2 = p.twiny()
#ay2.xaxis.tick_top()
#ay2.yaxis.tick_right()
#blmin, blmax = 0, 0.13 / dk_du(C.pspec.f2z(fq))
#p.axis([blmin,blmax,0,.6])
#p.xlabel('Baseline Length (Wavelength)',fontsize=16)
#title = r'$P(k)\ [{\rm mK}^2 (h^{-1}{\rm Mpc})^{3}]$'
#p.text(.5,1.15,title,horizontalalignment='center',transform=ay2.transAxes)


#p.savefig('data-wedge.eps', format='eps', dpi=250)
#p.savefig('data-wedge-gray.eps', format='eps', dpi=500)
p.show()
