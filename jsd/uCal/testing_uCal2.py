import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import scipy
import scipy.constants
import aipy as a
import optparse, sys, os
import capo
from scipy.sparse import csr_matrix

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

dataFiles = ['./Data/zen.2456679.49577.xx.npz']
#dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith('xx.npz') and file.startswith('zen.2456679.')]
#dataFiles = ["./Data/" + file for file in os.listdir("./Data") if file.endswith('xx.uvcRRE.npz') and file.startswith('zen.2456943')]
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
bls = getBaselines(freqs, intSeps, dataFiles[0], calFile=calFile)#[2:]
bl2SepDict = {bl: np.asarray([sep,0.0]) for bl,sep in zip(bls,separations)}
redundancyDict = blRedundancyDictFromCalFile(calFile=calFile)
chans = range(len(freqs))

flags = loadFlags(dataFiles, flagSuffix, flagThreshold)
data, samples = loadVisibilitiesAndSamples(dataFiles, pol, bls, redundancyDict, flags)


#%%


#sepBLPairs = sorted([[np.asarray(bl2SepDict[bl]), bl] for bl in bls], key=lambda pair: np.linalg.norm(pair[0])) #sorts by length


#Need to build a dict that converts chan and bl to uBin number


#chans = np.nonzero(np.sum(data[bls[0]]==0,axis=0)==0)[0]

def renormalize(betas, Sigmas):
    betaAbsMean = np.mean(np.abs(betas))
    betaAngleMean = np.mean(np.angle(betas))
    return betas / betaAbsMean * np.exp(-1.0j * betaAngleMean), Sigmas * betaAbsMean / np.exp(-1.0j * betaAngleMean)

allBetas = np.zeros((len(data[bls[0]]),len(freqs)), dtype=complex)
allBetasLog = np.zeros((len(data[bls[0]]),len(freqs)), dtype=complex)
allSigmas = []
for t in range(1):
#for t in range(len(data[bls[0]])):    
    chans = np.nonzero(data[bls[0]][t,:])[0]
    for toRemove in [15, 32, 34, 33, 186, 184, 181, 14, 16, 17, 18, 19, 74, 102, 128, 130, 166, 168, 169, 178, 182]:
        try: 
            chans = np.delete(chans, np.where(chans==toRemove)[0][0])
        except:
            print "No channel", toRemove
    uBinDict = {(bl,chan): np.floor((freqs[chan]*bl2SepDict[bl] - freqs[chans[0]]*np.array([separations[0],0]))/uBinSize).astype(int)[0] for bl in bls for chan in chans}    
    nuBins = np.max(np.unique(uBinDict.values()))+1
    nVis = len(bls)*len(chans)
    n = 0 
    Acoeffs, rowIndices, colIndices = [np.zeros(nVis*2) for i in range(3)]
    for i,chan in enumerate(chans):
        for bl in bls:
            rowIndices[2*n:2*n+2] = n
            colIndices[2*n:2*n+2] = [i, len(chans) + uBinDict[(bl,chan)]]
            Acoeffs[2*n:2*n+2] = [1.0, 1.0]
            n += 1
    
    logcalA = csr_matrix((Acoeffs,(rowIndices,colIndices)))
    logcalAtA = logcalA.conjugate().transpose().dot(logcalA).toarray()
    
    y = np.asarray([np.log(data[bl][t][chan]) for chan in chans for bl in bls])
    xhatReal = np.linalg.pinv(logcalAtA).dot((logcalA.conjugate().T).dot(np.real(y)))
    xhatImag = np.linalg.pinv(logcalAtA).dot((logcalA.conjugate().T).dot(np.imag(y)))
    betas = np.exp(xhatReal + 1.0j*xhatImag)[0:len(chans)]
    Sigmas = np.exp(xhatReal + 1.0j*xhatImag)[len(chans):]
    
    betas, Sigmas = renormalize(betas, Sigmas)
    allBetasLog[t, chans] = betas
    allSigmas.append(Sigmas)
#    continue
    
    
    
    noiseVariance = np.asarray([np.abs(beta)**2 / redundancyDict[bl]  for beta in betas for bl in bls])
    errors = np.asarray([betas[i] * Sigmas[uBinDict[bl,chans[i]]] - data[bl][t][chans[i]] for i in range(len(chans)) for bl in bls])
    noise = csr_matrix((np.append(noiseVariance**-1,noiseVariance**-1), (np.arange(2*nVis), np.arange(2*nVis))))
    chiSqPerDoF = np.average(np.append(np.real(errors)**2, np.imag(errors)**2) * noise)    
    print 'logcal', chiSqPerDoF
    
    for iteration in range(100):
#        noiseVariance = np.asarray([np.abs(beta)**2 / redundancyDict[bl]  for beta in betas for bl in bls])
        n = 0 
        coeffs, rowIndices, colIndices = [np.zeros(nVis*8) for i in range(3)]
        for i,chan in enumerate(chans):
            for bl in bls:
                beta0, Sigma0 = betas[i], Sigmas[uBinDict[(bl,chan)]]
                rowIndices[8*n:8*n+4] = n
                rowIndices[8*n+4:8*n+8] = n + nVis
                colIndices[8*n:8*n+8] = 2*[i, len(chans) + uBinDict[(bl,chan)], nuBins+len(chans) + i, nuBins + 2*len(chans) + uBinDict[(bl,chan)]]
                coeffs[8*n:8*n+8] = [np.real(beta0*Sigma0), np.real(beta0), -np.imag(beta0*Sigma0), -np.imag(beta0), 
                                        np.imag(beta0*Sigma0), np.imag(beta0), np.real(beta0*Sigma0), np.real(beta0)]
                n += 1
        A = csr_matrix((coeffs,(rowIndices,colIndices)))
        Ninv = csr_matrix((np.append((noiseVariance)**-1,(noiseVariance)**-1), (np.arange(2*nVis), np.arange(2*nVis))))
        AtNinvA = (A.conjugate().transpose().dot(Ninv)).dot(A).toarray()
        
        deltas = np.asarray([-betas[i] * Sigmas[uBinDict[bl,chan]] + data[bl][t][chan] for i,chan in enumerate(chans) for bl in bls])
        xHat = np.linalg.pinv(AtNinvA).dot(A.T.conjugate().dot(Ninv.dot(np.concatenate((np.real(deltas),np.imag(deltas))))))
        alpha = .5
        newBetas = np.asarray([betas[i] * (1.0 + alpha*(xHat[i] + 1.0j*xHat[i+nuBins+len(chans)])) for i in range(len(chans))])
        newSigmas = np.asarray([Sigmas[i] + alpha*(xHat[len(chans) + i] + 1.0j*xHat[i+nuBins+2*len(chans)]) for i in range(nuBins)])
        relativeChange = np.linalg.norm(newBetas - betas)/np.linalg.norm(betas)
        betas, Sigmas = newBetas, newSigmas
        betas, Sigmas = renormalize(betas, Sigmas)
        errors = np.asarray([betas[i] * Sigmas[uBinDict[bl,chans[i]]] - data[bl][t][chans[i]] for i in range(len(chans)) for bl in bls])
        noise = csr_matrix((np.append(noiseVariance**-1,noiseVariance**-1), (np.arange(2*nVis), np.arange(2*nVis))))
        chiSqPerDoF = np.average(np.append(np.real(errors)**2, np.imag(errors)**2) * noise)
        print iteration, chiSqPerDoF
        if relativeChange < 1e-3: break

        
        #%%
#        errorsByChan = {chan: [] for chan in chans}
#        n = 0
#        for i,chan in enumerate(chans):
#            for bl in bls:
#                errorsByChan[chan].append(errors[n]**2/noiseVariance[n])
#                n += 1
#        chanErrors = [np.mean(errorsByChan[chan])**.5 for chan in chans]
#        plt.figure(3); plt.clf()
#        plt.plot(chans, chanErrors,'.')
#            
        
        #%%
#        def lincalScatter(x,y,color=None, figNum=100, xs='log', ys='log', title='', clear=True, xl='', yl=''):
#            plt.figure(figNum); 
#            if clear: plt.clf()
#            if color is not None: plt.scatter(x,y,c=color)
#            else: plt.scatter(x,y)
#            plt.yscale(ys); plt.xscale(xs)
#            plt.ylim([.9*np.min(y), 1.1*np.max(y)]); plt.xlim([.9*np.min(x), 1.1*np.max(x)])
#            plt.xlabel(xl); plt.ylabel(yl); plt.title(title)
#        thisData = [data[bl][t][chans[i]] for i in range(len(chans)) for bl in bls]
#        thisModel = [betas[i] * Sigmas[uBinDict[bl,chans[i]]] for i in range(len(chans)) for bl in bls]
#        thisubin = [uBinDict[(bl,chan)] for i in range(len(chans)) for bl in bls]
#        thischan = [i for i in range(len(chans)) for bl in bls]
#        lincalScatter(np.abs(errors), np.abs(noiseVariance)**.5, figNum=2, color = thischan)
#        plt.plot([0,1],[0,1],'k--')
    
    allBetas[t, chans] = betas

#%%
plt.figure(1); plt.clf()
#plt.plot(np.mean(np.abs(allBetas.T),axis=1))
plt.plot(np.abs(allBetas.T),'.-')


plt.figure(11); plt.clf()
#plt.plot(np.mean(np.abs(allBetasLog.T),axis=1))
plt.plot(np.abs(allBetasLog.T),'.-')
#%%



#%%
"""
Psuedocode:
	For each time:
		perform logcal
		develop a noise variance
		perform lincal
"""

#TODO: consider discretizing u in terms of cos(u) modes rather than discrete points

#%%
plt.figure(100); plt.clf()
plt.plot(np.abs(data[bls[0]]).T)
