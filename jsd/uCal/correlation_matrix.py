import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import scipy
import scipy.constants
import aipy as a
import optparse, sys, os
import capo
import glob


def getBaselines(freqs, intSeps, exampleDataFile, calFile='psa6622_v003'):
    aa = a.cal.get_aa(calFile, freqs)
    blstr, _, bl2sep = zsa.grid2ij(aa.ant_layout)
    m, _, mdl,_= omni.from_npz(exampleDataFile)
    bls = [(0,i) for i in intSeps]
    for ij in mdl['xx'].keys():
        bl = a.miriad.ij2bl(ij[0],ij[1]) 
        try:
            s = tuple(map(int,bl2sep[bl].split(',')))
        except(KeyError):
            continue
        if s in bls and s[0] == 0:
            bls[s[-1]-1] = ij
    return bls

def blRedundancyDictFromCalFile(calFile = 'psa6622_v003', verbose = False):
    """ Uses Omnical to figure out the number of redundant baselines for a given separation (as demarkated by a representative baseline tuple). 
    Returns redundancyDict and also saves it as an instance variable."""
    aa = a.cal.get_aa(calFile,np.array([.15]))
    if verbose: print 'Now generating baseline redundancy dictionary from calfile...'
    ex_ants = [5,7,8,15,16,17,24,26,27,28,29,37,38,46,48,50,51,53,55,63,68,69,72,74,76,77,82,83,84,85,92,107,110]
    info = capo.omni.aa_to_info(aa, pols=['x'], ex_ants=ex_ants)
    reds = info.get_reds()
    redundancyDict = {}
    #This is overkill...it maps all baselines to their redundancy, not just the ones we need. But that's fine.
    for red in reds:
        for bl in red: redundancyDict[bl] = len(red)
    return redundancyDict

def loadFlags(dataFiles, flagSuffix, flagThreshold):
    """This function loads flags from xrfi and flags channels completely that are flagged more than the threshold fraction of the time"""
    if flagSuffix is None: return None
    flags = []
    for fl in dataFiles:
        flag = np.load(fl.replace('.npz', flagSuffix))
        for t in range(len(flag.keys())-1): flags.append(flag[str(t)])
    flags = np.asarray(flags)
    for chan in range(flags.shape[1]): 
        if np.sum(flags[:,chan]) > flagThreshold * flags.shape[0]: flags[:,chan] = 1
    return flags

def loadVisibilitiesAndSamples(dataFiles, pol, blList, redundancyDict, flags):
    """Loads visibilities and samples from omnical npzs. Applies flags."""
    print 'Now reading all data files...'
    data = {}
    samples = {}
    jds = []
    lsts = []
    files = []
    conjugate = [(0,101), (0,62), (0,100), (0,97), (12,43), (57,64)]
    for fl in dataFiles:
        print '    Reading %s'%fl
        meta,_,mdl,_ = omni.from_npz(fl)
        jd = meta['jds']
        lst = meta['lsts']
        jds.append(jd)
        lsts.append(lst)
        files.append(fl.split('/')[-1]) #ignore the folder
        for b in bls:
            if b in conjugate: _d = np.conj(mdl[pol][b])
            else: _d = mdl[pol][b]
            if b not in data.keys(): data[b] = _d
            else: data[b] = np.concatenate((data[b],_d), axis=0)
            samples[b] = redundancyDict[b] * np.ones(data[b].shape)
            samples[b][data[b]==0] = 0
    if flags is not None:
        for b in bls: data[b][flags==1] = 0 #apply flags
        for b in bls: samples[b][flags==1] = 0 #apply flags
    return data, samples

dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith('xx.npz') and file.startswith('zen.2456679.')]
#dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith('xx.npz') and file.startswith('zen.2456679.4')]
#dataFiles = ['./Data/zen.2456679.49577.xx.npz']
nBootstraps = 8
verbose = True
pol = 'xx'
flagSuffix, flagThreshold = None, .25
spectralIndex = 2.55
freqs = np.arange(.1, .2, .1/203)
chans = range(len(freqs))
calFile = 'psa6622_v003'
uMin, uMax, maxDeltau, uBinSize, duBinSize = 25.0, 125, .5, .5, 5.0/203
manualChannelFlags = [61, 62, 74, 178, 166, 169, 168, 102]#[14, 16, 17, 18, 19, 74, 102, 128, 130, 166, 168, 169, 178, 182, 184, 186]#
loadPickle = False
outfilePrefix = '.'

#PAPER 128 Setup
intSeps = np.arange(1,16) #* 15m  = baseline lengths
dx = 15.0 / scipy.constants.c * 1e9
chan2FreqDict = {chan: freqs[chan] for chan in chans}
separations = dx*intSeps
if calFile is None: calFile = 'psa6622_v003'
bls = getBaselines(freqs, intSeps, dataFiles[0], calFile=calFile)
bl2SepDict = {bl: np.asarray([sep,0.0]) for bl,sep in zip(bls,separations)}
redundancyDict = blRedundancyDictFromCalFile(calFile=calFile)
nChan = len(freqs)

#%%
flags = loadFlags(dataFiles, flagSuffix, flagThreshold)
data, samples = loadVisibilitiesAndSamples(dataFiles, pol, bls, redundancyDict, flags)

#%%
plt.figure(1); plt.clf()
blKeys = [item[1] for item in sorted({value[0]: key for key,value in bl2SepDict.items()}.items())]
allCorrs = np.zeros((nChan*len(blKeys),nChan*len(blKeys)))
for i in range(len(blKeys)):
    for j in range(len(blKeys)):
#        plt.subplot(len(blKeys), len(blKeys), i*len(blKeys)+j+1)
        allCorrs[i*nChan:(i+1)*nChan, j*nChan:(j+1)*nChan] = np.corrcoef(data[blKeys[i]].T, data[blKeys[j]].T)[nChan:,0:nChan]
plt.imshow(np.abs(allCorrs), interpolation='none')
#plt.imshow(np.abs(np.corrcoef(data[blKeys[i]].T, data[blKeys[j]].T)[nChan:,0:nChan]),interpolation='none')
plt.colorbar()