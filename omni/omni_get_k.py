#! /usr/bin/env python
import numpy as n
import aipy as a
import optparse,sys
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import pylab as p, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])


pol='xx'
freqs = n.arange(.1,.2,.1/203)
aa = a.cal.get_aa('psa6622_v003', freqs)
blstr, _, bl2sep = zsa.grid2ij(aa.ant_layout)
m, _, mdl,_= omni.from_npz(args[0])
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
        bls[s[-1]-1] = ij

#print bls
conjugate = [(0,101), (0,62), (0,100), (0,97), (12,43), (57,64)]

#kmodes = n.einsum('j,i', freqs, seps) #these are all the kmodes measured. Need to align them.

fdict = [{} for i in range(800)] #for a frequency, the seps and frequency that are redundant
du = 15*15*(freqs[2]-freqs[1])/3e8 * 1e9 #u-range covered b/w neighboring freqs
us = (n.arange(400)+.5)*du #.5*du*N 
errors = n.ones(shape=(len(freqs),len(freqs)))
cut = 60/15.*du
print 'Finding baseline/frequency pairs'
for ch1,fq1 in enumerate(freqs):
    #print ch1,
    for ch2,fq2 in zip(range(ch1+1,len(freqs)),freqs[ch1+1:]):
        if ch1 == ch2: continue
        err = []
        blp = []
        bestsofar = [0,(0,0),(0,0),0,100000.] #u-bin, freqpair, blpair, u-value,u-error
        for s1 in n.arange(15,0,-1):
            for s2 in n.arange(s1-1,0,-1):
                if s1==s2: continue
                if n.abs(fq1*s1 - fq2*s2)*15*1e9/3e8 < min(cut,bestsofar[-1]):
                    u = (fq1*s1 + fq2*s2)*15/2.*1e9/3e8
                    #fdict[int(n.floor(u/du))][(ch1,ch2)] = [(bls[s1-1],bls[s2-1]),  n.abs(fq1*s1 - fq2*s2)*15*1e9/3e8]
                    bestsofar[0] = int(n.floor(u/du))
                    bestsofar[1] = (ch1,ch2)
                    bestsofar[2] = (bls[s1-1],bls[s2-1])
                    bestsofar[3] = u
                    bestsofar[4] = n.abs(fq1*s1 - fq2*s2)*15*1e9/3e8
        fdict[bestsofar[0]][bestsofar[1]] = [bestsofar[2],bestsofar[3],bestsofar[4]]

for i,d in enumerate(fdict):
    if len(d) <= 1 : fdict[i] = {}
#                err.append(n.abs(fq1*s1 - fq2*s2))
#                blp.append((s1,s2))
#        errors[ch1,ch2] = n.min(err)
#        if errors[ch1,ch2]<cut:
#            fdict[(ch1,ch2)] = n.asarray(blp)[n.where(err == n.min(err))]
#        else: continue
#print 

#print n.sum(errors - errors.T)
p.imshow(errors, origin='lower', aspect='equal', extent=(freqs[0],freqs[-1],freqs[0],freqs[-1]))

#read in data
data = {}; 
jds = []
lsts = []
files = []
prefix = '5' #change if don't want to overwrite filenames
for fl in args:
    if os.path.exists(args[-1][:-3]+'ucal' + prefix + '.npz'):
        print '   '+args[-1][:-3]+'ucal' + prefix + '.npz exists. Skipping...'
        exit()
    print 'Reading %s'%fl
    meta,_,mdl,_ = omni.from_npz(fl)
    jd = meta['jds']
    lst = meta['lsts']
    jds.append(jd)
    lsts.append(lst)
    files.append(fl)
    for b in bls:
        if b in conjugate:
            _d = n.conj(mdl[pol][b])
        else:
            _d = mdl[pol][b]
        if b not in data.keys(): data[b] = _d
        else: 
            data[b] = n.concatenate((data[b],_d), axis=0)

reddata = {}
reddata['files'] = files
for i,d in enumerate(fdict):
    for ch1,ch2 in d.keys():
        b1,b2 = d[(ch1,ch2)][0]
        #getting rid of partially flagged channels
        w1 = data[b1][:,ch1] != 0
        w2 = data[b2][:,ch2] != 0
        w = n.logical_and(w1,w2)
        if n.all(n.logical_not(w)): continue
        avis = n.average(data[b1][:,ch1]*n.conj(data[b2][:,ch2]), weights=w)
#        if not avis: continue
        reddata[str((ch1,ch2))] = [ (b1,b2), i, avis ]
print '   Saving file ' + fl[:-3] + 'ucal' + prefix + '.npz' #name is last file
n.savez(fl[:-3] + 'ucal' + prefix + '.npz',**reddata)

#print len(reddata)

#import IPython
#IPython.embed()

x,y = n.meshgrid(freqs, seps)
z = x*y
cs = p.contour(x,y,z, 100)
p.scatter(x,y,s=1.5)
p.grid(True)
#p.show()
