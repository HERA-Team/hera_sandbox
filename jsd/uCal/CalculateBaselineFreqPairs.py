#! /usr/bin/env python

#Based on Zaki and Carina's omni_get_k.py
import numpy as np
import aipy as a
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import optparse, sys, os
import cPickle as pickle

pol='xx'
seps = np.arange(1,16) #* 15m  = baseline lengths
dx = 15.0
separations = dx*seps
freqs = np.arange(.1,.2,.1/203)
uTolerance = 15.0 * 15 * (freqs[1]-freqs[0]) / 3e8 * 1e9  

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,dataFiles = o.parse_args(sys.argv[1:])

if len(sys.argv[1:]) == 0: #This option is for testing
    #dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("uvcRRE.npz")] 
    dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith("2456943.57058.xx.uvcRRE.npz")] 
    
aa = a.cal.get_aa('psa6622_v003', freqs)
blstr, _, bl2sep = zsa.grid2ij(aa.ant_layout)
m, _, mdl,_= omni.from_npz(dataFiles[0])

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


baselineFreqPairs = {}
print 'Now finding all baseline/frequency pairs...'
for ch1,f1 in enumerate(freqs):
    for ch2,f2 in zip(range(ch1+1,len(freqs)),freqs[ch1+1:]):
        for blIndex2,sep2,pair2 in zip(range(len(separations)), separations, bls):
            for blIndex1,sep1,pair1 in zip(range(blIndex2+1,len(separations)), separations[blIndex2+1:], bls[blIndex2+1:]):
                if abs(f1*sep1 - f2*sep2) * 1.0e9 / 3.0e8 < uTolerance:
                    u = (f1*sep1 + f2*sep2)/2.0*1e9/3e8
                    if not baselineFreqPairs.has_key((f1, f2, u)) and not baselineFreqPairs.has_key((f2, f1, u)):
                        baselineFreqPairs[(f1, f2, u)] = {'separations': (sep1, sep2), 'du': abs(f1*sep1 - f2*sep2)*1.0e9/3.0e8, 'chans': (ch1,ch2), 'blpairs': (pair1,pair2)}
print "   " + str(len(baselineFreqPairs)) + " baseline/frequency pairs identified with delta u < " + str(uTolerance)

# DATA LOADING STUFF

print 'Now reading all data files...'
#read in data
data = {}; 
jds = []
lsts = []
files = []
prefix = '' #change if don't want to overwrite filenames
conjugate = [(0,101), (0,62), (0,100), (0,97), (12,43), (57,64)]
for fl in dataFiles:
    print '   Reading %s'%fl
    meta,_,mdl,_ = omni.from_npz(fl)
    jd = meta['jds']
    lst = meta['lsts']
    jds.append(jd)
    lsts.append(lst)
    files.append(fl.split('/')[-1]) #ignore the folder
    for b in bls:
        if b in conjugate:
            _d = np.conj(mdl[pol][b])
        else:
            _d = mdl[pol][b]
        if b not in data.keys(): data[b] = _d
        else: 
            data[b] = np.concatenate((data[b],_d), axis=0)


#%%
# #SAVE RESULTS
if 'files' in baselineFreqPairs: del baselineFreqPairs['files']
for (f1,f2,u), entry in baselineFreqPairs.items():
    b1, b2 = entry['blpairs']
    ch1, ch2 = entry['chans']
    w1 = data[b1][:,ch1] != 0 #getting rid of partially flagged channels
    w2 = data[b2][:,ch2] != 0
    w = np.logical_and(w1,w2)
    if np.all(np.logical_not(w)): del baselineFreqPairs[(f1,f2,u)]
    else: 
        visProds = (data[b1][:,ch1]*np.conj(data[b2][:,ch2]))[w]
        baselineFreqPairs[(f1,f2,u)]['avgvis'] = np.average(visProds)
        baselineFreqPairs[(f1,f2,u)]['nIntegrations'] = len(visProds)
        baselineFreqPairs[(f1,f2,u)]['averageDeltaProd'] = np.average(np.abs(visProds[0:-2] - visProds[1:-1]))
        baselineFreqPairs[(f1,f2,u)]['stdDeltaProd'] = np.std(visProds[0:-2] - visProds[1:-1])

baselineFreqPairs['files'] = files
print 'Now saving file ' + dataFiles[-1][:-3] + 'ucal' + prefix + '.p with ' + str(len(baselineFreqPairs)-1) + ' entries' #name is last file
pickle.dump(baselineFreqPairs, open(dataFiles[-1][:-3] + 'ucal' + prefix + '.p', 'wb'))

