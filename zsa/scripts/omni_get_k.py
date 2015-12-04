#! /usr/bin/env python -u
import numpy as n
import aipy as a
import optparse,sys
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import pylab as p

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])


pol='xx'
freqs = n.arange(.1,.2,.1/203)
aa = a.cal.get_aa('psa6622_v003', freqs)
blstr, _, bl2sep = zsa.grid2ij(aa.ant_layout)
_, _, mdl, _, jds, lsts, freqs = omni.from_npz(args[0])
seps = n.arange(1,16) #* 15m  = baseline lengths
bls = [(0,i) for i in seps]
sep2ij = {}
for ij in mdl['xx'].keys():
    bl = a.miriad.ij2bl(ij[0],ij[1]) 
    try:
        s = tuple(map(int,bl2sep[bl].split(',')))
    except(KeyError):
        continue
    if s in bls and s[0] == 0:
#        print bl, ij, s
        bls[s[-1]-1] = ij

print bls
conjugate = [(0,101), (0,62), (0,100), (0,97), (12,43), (57,64)]

kmodes = n.einsum('j,i', freqs, seps) #these are all the kmodes measured. Need to align them.

fdict = {} #for a frequency, the seps and frequency that are redundant.
errors = n.ones(shape=(len(freqs),len(freqs)))
thresholdcut = 18.0 #degree error. 18 =.01
cut = thresholdcut/360/1e9*3e8/15
print cut
for ch1,fq1 in enumerate(freqs):
    print ch1,
    for ch2,fq2 in enumerate(freqs):
#        if (ch1,ch2) not in fdict.keys(): fdict[(ch1,ch2)] = None
        if ch1 == ch2: continue
        err = []
        blp = []
        for s1 in n.arange(1,16):
            for s2 in n.arange(1,16):
#                if s1==s2: continue
                err.append(n.abs(fq1*s1 - fq2*s2))
                blp.append((s1,s2))
        errors[ch1,ch2] = n.min(err)
        if errors[ch1,ch2]<cut:
            fdict[(ch1,ch2)] = n.asarray(blp)[n.where(err == n.min(err))]
        else: continue

print 'kjdsfakjsnfa'
print 
print n.sum(errors - errors.T)
#arp.waterfall(errors,mode='lin')
p.imshow(errors, origin='lower', aspect='equal', extent=(freqs[0],freqs[-1],freqs[0],freqs[-1]))
#ns, ds = p.hist(n.flatten(errors), bins=100)
#print ns.shape, ds.shape
#p.colorbar(shrink=.5)

#read in data
data = {}; 
jds = []
lsts = []
for fl in args:
    print 'reading %s'%fl
    _,_,mdl,_,jd,lst,fq = omni.from_npz(fl)
    jds.append(jd)
    lsts.append(lst)
    for b in bls:
        if b in conjugate:
#            print b
            _d = n.conj(mdl[pol][b])
        else:
            _d = mdl[pol][b]
        if b not in data.keys(): data[b] = _d
        else: 
            data[b] = n.concatenate((data[b],_d), axis=0)


reddata = {}
for ch1,ch2 in fdict.keys():
    try:
        ss = fdict[(ch1,ch2)][0]
    except(TypeError): 
        continue
    s1,s2 = ss[0],ss[1]
    b1,b2 = bls[s1-1], bls[s2-1]
    avis = n.mean(data[b1][:,ch1]*n.conj(data[b2][:,ch2]))
    if not avis: continue
    else:
        reddata[str((ch1,ch2))] = [ (b1,b2), avis ]


n.savez('ucal.npz',**reddata)
                
    

        
#arp.waterfall(n.array(data[b]), mode='lin')
#p.show()
    


x,y = n.meshgrid(freqs, seps)
#print y
z = x*y
#print z
cs = p.contour(x,y,z, 100)
#print cs
p.scatter(x,y,s=1.5)
#p.colorbar(shrink=.5)
#p.clabel(cs, inline=1)
p.grid(True)
p.show()
