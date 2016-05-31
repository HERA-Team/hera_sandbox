#import healpy as hp, numpy as n

#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p

#@p.ion()
#fqs = n.linspace(.1,.2,203)
fq = .15
bl1, bl2 = (0,26),(0,26)
N = 5   #number of universes to average over
aa = a.cal.get_aa('psa6240_v003', n.array([fq]))
h = a.healpix.HealpixMap(nside=64)
lmax = 300
mmax = lmax
ALM,CL = {},n.zeros(lmax)
for i in xrange(N):
    print i
    sky = n.random.normal(size=h.map.size)
    h.map = sky # assume sky is in eq coord
    alm = h.to_alm(lmax,mmax).get_data()
    ind = 0
    for l in n.arange(0,lmax+1):
    	cl=0
    	for m in range(0,l+1):
    		if m>=1: cl+=2*(n.conj(alm[ind])*alm[ind])/(2*l+1)   #<almal'm'>=delll'delmm'Cl
    		else: cl+=(n.conj(alm[ind])*alm[ind])/(2*l+1)
    	ALM[l] = ALM.get(l,[])+[cl]
    for L in ALM.keys(): 
    	#print L, n.mean(ALM[L]).real
    	CL[int(L-1)] += n.mean(ALM[L]).real/N
    print CL


#fig = p.figure(1)
#hp.mollview(h, min=-1, max=1, title='Unmasked map', fig=1, unit=r'$\Delta$T (mK)')
#p.show()

