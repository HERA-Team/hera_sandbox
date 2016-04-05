#! /usr/bin/env python

import numpy as np
import cPickle as pickle
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import scipy
import omnical
import time
import aipy
import capo
plt.close('all')

def SortedEigensystem(matrix):
    """Returns the eigensystem of the input matrix where eigenvalues and eigenvectors are sorted by descending absolute value."""
    evals,evecs = np.linalg.eig(matrix)
    indices = np.argsort(np.abs(evals))[::-1]   
    return evals[indices], evecs[:,indices]

#baselineFreqPairs = pickle.load(open('./Data/zen.2456943.65409.xx.uvcRRE.ucal.p','rb'))
baselineFreqPairs = pickle.load(open('./Data/zen.2456943.57058.xx.uvcRRE.ucal.p','rb'))

files = baselineFreqPairs['files']; del baselineFreqPairs['files']
uMinCutoff = 5
uMaxCutoff = 125
duThreshold = .25

#flaggedChannels = np.asarray([])
flaggedChannels = np.asarray([14, 55, 101, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186])
#TODO: find an automatic way to flag channels



for key in baselineFreqPairs.keys(): 
    if key[2] > uMaxCutoff: del baselineFreqPairs[key]
    elif baselineFreqPairs[key]['du'] > duThreshold: del baselineFreqPairs[key]
    elif key[2] < uMinCutoff: del baselineFreqPairs[key]
    elif (baselineFreqPairs[key]['chans'][0] == flaggedChannels).any(): del baselineFreqPairs[key]
    elif (baselineFreqPairs[key]['chans'][1] == flaggedChannels).any(): del baselineFreqPairs[key]
nPairs = len(baselineFreqPairs)
print 'There are ' + str(nPairs) + ' baseline frequency pairs after our cuts.'
allFreqs = np.unique([[key[0],key[1]] for key in baselineFreqPairs.keys()])
nFreqs = len(allFreqs)
allus = np.unique([key[2] for key in baselineFreqPairs.keys()])
alldus = np.unique([entry['du'] for key,entry in baselineFreqPairs.items()])

#Assign uBins for logcal
nuBins = 50
uMin, uMax = np.min(allus), np.max(allus)
deltau = (uMax-uMin)/(nuBins-1e-9)
for key in baselineFreqPairs.keys(): baselineFreqPairs[key]['uBin'] = int(np.floor((key[2]-uMin)/deltau))
alluBins, uBinWeights = np.unique([entry['uBin'] for key,entry in baselineFreqPairs.items()],return_counts=True)

uBinCenters = {n:[] for n in range(nuBins)}
for key,entry in baselineFreqPairs.items(): uBinCenters[entry['uBin']].append(key[2])
for uBin,uValues in uBinCenters.items(): uBinCenters[uBin] = np.average(uValues)
for uBin,uValues in uBinCenters.items(): 
    if np.isnan(uBinCenters[uBin]): del uBinCenters[uBin]
    

#%%    
#Assign duBins for logcal
duMin, duMax = np.min(alldus), np.max(alldus)
nduBins = int(np.round((duMax-duMin) / (.1/203/3e8*1e9*15))) + 1
deltadu = (duMax-duMin)/(nduBins-1+1e-9)
for key in baselineFreqPairs.keys(): baselineFreqPairs[key]['duBin'] = int(np.floor((baselineFreqPairs[key]['du']-duMin)/deltadu))
allduBins, duBinWeights = np.unique([entry['duBin'] for key,entry in baselineFreqPairs.items()],return_counts=True)

duBinCenters = {n:[] for n in range(nduBins)}
for key,entry in baselineFreqPairs.items(): duBinCenters[entry['duBin']].append(entry['du'])
for duBin,duValues in duBinCenters.items(): duBinCenters[duBin] = np.average(duValues)
for duBin,duValues in duBinCenters.items(): 
    if np.isnan(duBinCenters[duBin]): del duBinCenters[duBin]

#%% Perform logcal (assuming N = I)

print '\nNow performing Logcal...'

freq2Col = {freq: col for col,freq in zip(range(0,nFreqs),allFreqs)}
uBin2Col = {uBin: col for col,uBin in zip(range(nFreqs,nFreqs+nuBins),alluBins)}
duBin2Col = {duBin: col for col,duBin in zip(range(nFreqs+nuBins,nFreqs+nuBins+nduBins),allduBins)}
col2freq = {col: freq for freq,col in freq2Col.items()}
col2uBin = {col: uBin for uBin,col in uBin2Col.items()}
col2duBin = {col: duBin for duBin,col in duBin2Col.items()}

Acoeffs, Bcoeffs, rowIndices, colIndices = np.zeros(nPairs*4), np.zeros(nPairs*4), np.zeros(nPairs*4), np.zeros(nPairs*4)
for n,((f1,f2,u),entry) in enumerate(baselineFreqPairs.items()):
    rowIndices[4*n:4*n+4] = n
    colIndices[4*n:4*n+4] = [freq2Col[f1], freq2Col[f2], uBin2Col[entry['uBin']], duBin2Col[entry['duBin']]]
    Acoeffs[4*n:4*n+4] = [1.0, 1.0, 1.0, 1.0]
    Bcoeffs[4*n:4*n+4] = [1.0, -1.0, 1.0, 1.0]

Amatrix = csr_matrix((Acoeffs,(rowIndices,colIndices)),shape=(nPairs+2, nFreqs+nuBins+nduBins))
Bmatrix = csr_matrix((Bcoeffs,(rowIndices,colIndices)),shape=(nPairs+2, nFreqs+nuBins+nduBins))
for freq in allFreqs:
    Amatrix[nPairs+0, freq2Col[freq]] = 1.0/nFreqs #Real parts average to 0 
    Bmatrix[nPairs+0, freq2Col[freq]] = 1 #Imag parts add to 0
for duBin in allduBins:
    Amatrix[nPairs+1, duBin2Col[duBin]] = 1.0/nduBins #Real parts average to 0
    Bmatrix[nPairs+1, duBin2Col[duBin]] = 1 #Imag parts add to 0
    

AtAlog = Amatrix.conjugate().transpose().dot(Amatrix).toarray()
BtBlog = Bmatrix.conjugate().transpose().dot(Bmatrix).toarray()

print 'Zero Eigenvalues: ', [np.sum(np.abs(np.linalg.eigvals(XtransX)<1e-10)) for XtransX in [AtAlog, BtBlog]]

yReal = np.append(np.real(np.log(np.asarray([item[1]['avgvis'] for item in baselineFreqPairs.items()]))),[0,0])
xhatReal = np.linalg.pinv(AtAlog).dot(Amatrix.conjugate().T.dot(yReal))
yImag = np.append(np.imag(np.log(np.asarray([item[1]['avgvis'] for item in baselineFreqPairs.items()]))),[0.0,0.0])
xhatImag = np.linalg.pinv(BtBlog).dot(Bmatrix.conjugate().T.dot(yImag))
logcalBandpass = np.exp(xhatReal + 1.0j*xhatImag)


plt.figure(1); plt.clf()
plt.plot(allFreqs,[np.abs(logcalBandpass[freq2Col[freq]]) for freq in allFreqs],'.')
plt.title('Logcal abs(beta)')
plt.xlabel('f (GHz)')

plt.figure(2); plt.clf()
plt.plot([uBinCenters[uBin] for uBin in alluBins],[np.abs(logcalBandpass[uBin2Col[uBin]]) for uBin in alluBins],'.')
plt.title('Logcal abs(Sigma(u))')
plt.xlabel('u')



#%% Compute Noise Covariance from Omnical

#Hardcoded for now: array information hard coded for now!
aa = aipy.cal.get_aa('psa6622_v003',np.array([.15]))
print 'Getting reds from calfile'
ex_ants = [5,7,8,15,16,17,24,26,27,28,29,37,38,46,48,50,51,53,55,63,68,69,72,74,76,77,82,83,84,85,92,107,110]
info = capo.omni.aa_to_info(aa, pols=['x'], ex_ants=ex_ants)
reds = info.get_reds()
redundancyDict = {}
for red in reds:
    for bl in red: redundancyDict[bl] = len(red)

spectralIndex = -2.55
freq2Col = {freq: col for col,freq in zip(range(0,nFreqs),allFreqs)}
#noiseCovDiag = np.asarray([(f1*f2)**(2*spectralIndex) / (1.0 * redundancyDict[entry['blpairs'][0]] * redundancyDict[entry['blpairs'][1]] * entry['nIntegrations']) for ((f1,f2,u),entry) in baselineFreqPairs.items()])
noiseCovDiag = np.asarray([(np.abs(logcalBandpass[freq2Col[f1]])**2 * np.abs(logcalBandpass[freq2Col[f2]])**2) / (1.0 * redundancyDict[entry['blpairs'][0]] * redundancyDict[entry['blpairs'][1]] * entry['nIntegrations']) for ((f1,f2,u),entry) in baselineFreqPairs.items()])
visCorrs = np.asarray([entry['avgvis'] for key,entry in baselineFreqPairs.items()])
noiseCovDiag *= np.var(visCorrs) / np.mean(noiseCovDiag) * 8e-6 * 750.0
#TODO: need to renormalize the noise


#%% Normal Lincal

print '\nNow performing Lincal with Discrete Sigma(u) and D(du)...'

betas = {freq: logcalBandpass[i] for i,freq in enumerate(allFreqs)}
uBin2Col = {uBin: col for col,uBin in zip(range(nFreqs,nFreqs+nuBins),alluBins)}
Sigmas = {uBin: logcalBandpass[uBin2Col[uBin]] for uBin in alluBins}
duBin2Col = {duBin: col for col,duBin in zip(range(nFreqs+nuBins,nFreqs+nuBins+nduBins),allduBins)}
Ds = {duBin: logcalBandpass[duBin2Col[duBin]] for duBin in allduBins}
visCorrs = np.asarray([entry['avgvis'] for key,entry in baselineFreqPairs.items()])

freq2Col = {freq: col for col,freq in zip(range(0,2*nFreqs,2),allFreqs)}
uBin2Col = {uBin: col for col,uBin in zip(range(2*nFreqs,2*nFreqs+2*nuBins,2),alluBins)}
duBin2Col = {duBin: col for col,duBin in zip(range(2*nFreqs+2*nuBins,2*nFreqs+2*nuBins+2*nduBins,2),allduBins)}
col2freq = {col: freq for freq,col in freq2Col.items()}
col2uBin = {col: uBin for uBin,col in uBin2Col.items()}
col2duBin = {col: duBin for duBin,col in duBin2Col.items()}

alpha = .3

for iteration in range(8):

    deltas = visCorrs - np.asarray([betas[f1]*np.conj(betas[f2])*Sigmas[entry['uBin']]*Ds[entry['duBin']] for ((f1,f2,u),entry) in baselineFreqPairs.items()])
    coeffs, rowIndices, colIndices = np.zeros(nPairs*16), np.zeros(nPairs*16), np.zeros(nPairs*16)        
    
    for n,((f1,f2,u),entry) in enumerate(baselineFreqPairs.items()):
        bbstarSD = betas[f1]*np.conj(betas[f2])*Sigmas[entry['uBin']]*Ds[entry['duBin']] #(beta)(beta^*)(Sigma)(D)
        bbstarS = betas[f1]*np.conj(betas[f2])*Sigmas[entry['uBin']] #(beta)(beta^*)(Sigma)
        bbstarD = betas[f1]*np.conj(betas[f2])*Ds[entry['duBin']] #(beta)(beta^*)(D)
        f1Col, f2Col, uCol, duCol = freq2Col[f1], freq2Col[f2], uBin2Col[entry['uBin']], duBin2Col[entry['duBin']]
        
        rowIndices[16*n:18*n+6] = n #the first 8 terms are on this row
        rowIndices[16*n+8:16*n+16] = n+nPairs #the next 8 terms are on the corresponding imaginary part row

        for i,colIndex in enumerate([f1Col,f2Col,f1Col+1,f2Col+1,uCol,uCol+1,duCol,duCol+1,f1Col,f2Col,f1Col+1,f2Col+1,uCol,uCol+1,duCol,duCol+1]):
            colIndices[16*n+i] = colIndex #these are eta1, eta2, phi1, phi2, sigma, psi, etc.
        
        coeffList = [np.real(bbstarSD), np.real(bbstarSD), -np.imag(bbstarSD), np.imag(bbstarSD), np.real(bbstarD), -np.imag(bbstarD), np.real(bbstarS), -np.imag(bbstarS),
                     np.imag(bbstarSD), np.imag(bbstarSD), np.real(bbstarSD), -np.real(bbstarSD), np.imag(bbstarD), np.real(bbstarD), np.imag(bbstarS), np.real(bbstarS)]
        for i,coeff in enumerate(coeffList): 
            coeffs[16*n+i] = coeff #these are the coefficients of those terms
    
    A = csr_matrix((coeffs,(rowIndices,colIndices)))
    #Ninv = csr_matrix((np.append(np.append(noiseCovDiag**-1,noiseCovDiag**-1),[1,1,1]), (np.arange(2*nPairs+3), np.arange(2*nPairs+3))))
    Ninv = csr_matrix((np.append(noiseCovDiag**-1,noiseCovDiag**-1), (np.arange(2*nPairs), np.arange(2*nPairs))))
    AtNinvA = (A.conjugate().transpose().dot(Ninv)).dot(A).toarray()
    if iteration == 0: print 'Zero Eigenvalues: ', len(AtNinvA) - np.linalg.matrix_rank(AtNinvA)

    deltasSplitReIm = np.append(np.real(deltas),np.imag(deltas))
    xHat = np.linalg.pinv(AtNinvA).dot(A.T.conjugate().dot(Ninv.dot(deltasSplitReIm)))

    newBetas = {freq: betas[freq]*(1+alpha*(xHat[freq2Col[freq]] + 1.0j*xHat[freq2Col[freq]+1])) for freq in allFreqs}
    newSigmas = {uBin: Sigmas[uBin] + alpha*(xHat[uBin2Col[uBin]] + 1.0j*xHat[uBin2Col[uBin]+1]) for uBin in alluBins}
    newDs = {duBin: Ds[duBin] + alpha*(xHat[duBin2Col[duBin]] + 1.0j*xHat[duBin2Col[duBin]+1]) for duBin in allduBins}
    betas = newBetas    
    Sigmas = newSigmas
    Ds = newDs

    #Renormalize    TODO: FIX NORMALIZATION AND PHASES IN THE LINEAR ALGEBRA ITSELF
    bandpassMean = np.average(np.asarray([np.abs(betas[f]) for f in allFreqs]))    
    for f in allFreqs: betas[f] = betas[f] / bandpassMean
    for key in Sigmas.keys(): Sigmas[key] = Sigmas[key] * bandpassMean**2
    bandpassAngleMean = np.average(np.asarray([np.angle(betas[f]) for f in allFreqs]))
    for f in allFreqs: betas[f] = betas[f] * np.exp(-1.0j*bandpassAngleMean)
    Dmean = np.average(np.asarray([np.abs(Ds[duBin]) for duBin in allduBins]))
    for key in Sigmas.keys(): Sigmas[key] = Sigmas[key] * Dmean
    for key in Ds.keys(): Ds[key] = Ds[key] / Dmean

    
    allChiSq = np.abs(visCorrs - np.asarray([betas[f1]*np.conj(betas[f2])*Sigmas[entry['uBin']]*Ds[entry['duBin']] for ((f1,f2,u),entry) in baselineFreqPairs.items()]))**2 / noiseCovDiag
    plt.figure(10); plt.clf()
    plt.plot(allFreqs,[np.abs(betas[freq]) for freq in allFreqs],'.')
    plt.title('Discrete Lincal abs(beta)')
    plt.xlabel('f (GHz)')
    plt.figure(11); plt.clf()
    plt.plot(uBinCenters.values(), [np.abs(Sigmas[uBin]) for uBin in alluBins],'.')
    plt.title('Discrete Lincal abs(Sigma(u))')
    plt.xlabel('u')
    plt.figure(12); plt.clf()
    plt.plot(allduBins, [np.abs(Ds[duBin]) for duBin in allduBins],'.')
    plt.title('Discrete Lincal abs(D(du))')
    plt.xlabel('du')
    plt.show() 
    plt.draw()
    plt.pause(0.05)
    print str(iteration) + ') mean chi^2 = ' + str(np.average(allChiSq))

binnedLincalAtNinvA = AtNinvA.copy()
binnedLincalBetas = betas.copy()
binnedLincalSigmas = Sigmas.copy()
binnedLincalDs = Ds.copy()
    
#%% Fit Bandpass to intialize lincal

realSigmaPolyOrder = 6
imagSigmaPolyOrder = realSigmaPolyOrder
uMean = np.mean(np.asarray([u for (f1,f2,u) in baselineFreqPairs.keys()]))
uBinCenterList = np.asarray([uBinCenters[key] for key in uBinCenters.keys()])
SigmaList = np.asarray([binnedLincalSigmas[key] for key in uBinCenters.keys()])
initialRealSigmaPolyCoeffs = np.polyfit(uBinCenterList/uMean, np.real(SigmaList), realSigmaPolyOrder, w=uBinWeights)
initialImagSigmaPolyCoeffs = np.polyfit(uBinCenterList/uMean, np.imag(SigmaList), imagSigmaPolyOrder, w=uBinWeights)

plt.figure(21); plt.clf()
plt.plot(uBinCenterList/uMean,np.real(SigmaList),'.')
plt.plot(uBinCenterList/uMean,np.polyval(initialRealSigmaPolyCoeffs,uBinCenterList/uMean))    
plt.plot(uBinCenterList/uMean,np.imag(SigmaList),'.')
plt.plot(uBinCenterList/uMean,np.polyval(initialImagSigmaPolyCoeffs,uBinCenterList/uMean))    

#%% Fit D(du) 

realDsPolyOrder = 4
imagDsPolyOrder = realDsPolyOrder
duMean = np.mean(np.asarray([entry['du'] for key,entry in baselineFreqPairs.items()]))
duBinCenterList = np.asarray([duBinCenters[key] for key in duBinCenters.keys()])
DsList = np.asarray([binnedLincalDs[key] for key in duBinCenters.keys()])
initialRealDsPolyCoeffs = np.polyfit(duBinCenterList/duMean, np.real(DsList), realDsPolyOrder, w=duBinWeights)
initialImagDsPolyCoeffs = np.polyfit(duBinCenterList/duMean, np.imag(DsList), imagDsPolyOrder, w=duBinWeights)

plt.figure(22); plt.clf()
plt.plot(duBinCenterList/duMean,np.real(DsList),'.')
plt.plot(duBinCenterList/duMean,np.polyval(initialRealDsPolyCoeffs,duBinCenterList/duMean))    
plt.plot(duBinCenterList/duMean,np.imag(DsList),'.')
plt.plot(duBinCenterList/duMean,np.polyval(initialImagDsPolyCoeffs,duBinCenterList/duMean))    


#%% Perform lincal with polynomial

print '\nNow performing Lincal with Polynomial Sigma(u) and D(du)...'

alpha = .5

betas = binnedLincalBetas
realSigmas = initialRealSigmaPolyCoeffs
imagSigmas = initialImagSigmaPolyCoeffs
realDs = initialRealDsPolyCoeffs
imagDs = initialImagDsPolyCoeffs
freq2Col = {freq: col for col,freq in zip(range(0,2*nFreqs,2),allFreqs)}
col2freq = {col: freq for freq,col in freq2Col.items()}


for iteration in range(8):
    betaBetaStarList = np.asarray([betas[f1]*np.conj(betas[f2]) for ((f1,f2,u),entry) in baselineFreqPairs.items()])    
    SigmasList = np.asarray([(np.polyval(realSigmas,u/uMean) + 1.0j*np.polyval(imagSigmas,u/uMean)) for ((f1,f2,u),entry) in baselineFreqPairs.items()])
    DsList = np.asarray([(np.polyval(realDs,entry['du']/duMean) + 1.0j*np.polyval(imagDs,entry['du']/duMean)) for ((f1,f2,u),entry) in baselineFreqPairs.items()])
    deltas = visCorrs - betaBetaStarList * SigmasList * DsList
    termsPerRow = 4*2 + (realSigmaPolyOrder+1)*2 + (imagSigmaPolyOrder+1)*2 + (realDsPolyOrder+1)*2 + (imagDsPolyOrder+1)*2
    coeffs, rowIndices, colIndices = np.zeros(nPairs*termsPerRow), np.zeros(nPairs*termsPerRow), np.zeros(nPairs*termsPerRow)
    realSigmaPolyCols = range(2*nFreqs,2*nFreqs+realSigmaPolyOrder+1)
    imagSigmaPolyCols = range(2*nFreqs+realSigmaPolyOrder+1,2*nFreqs+realSigmaPolyOrder+imagSigmaPolyOrder+2)
    realDsPolyCols = range(2*nFreqs+realSigmaPolyOrder+imagDsPolyOrder+2, 2*nFreqs+realSigmaPolyOrder+imagDsPolyOrder+realDsPolyOrder+3)
    imagDsPolyCols = range(2*nFreqs+realSigmaPolyOrder+imagDsPolyOrder+realDsPolyOrder+3, 2*nFreqs+realSigmaPolyOrder+imagDsPolyOrder+realDsPolyOrder+imagDsPolyOrder+4)
    for n,((f1,f2,u),entry) in enumerate(baselineFreqPairs.items()):
        du = entry['du']
        SigmaHere = (np.polyval(realSigmas,u/uMean) + 1.0j*np.polyval(imagSigmas,u/uMean))
        DHere = (np.polyval(realDs,du/duMean) + 1.0j*np.polyval(imagDs,du/duMean))
        bbstarSD = betas[f1]*np.conj(betas[f2])*SigmaHere*DHere #(beta)(beta^*)(Sigma)(D)
        bbstarS = betas[f1]*np.conj(betas[f2])*SigmaHere #(beta)(beta^*)(Sigma)
        bbstarD = betas[f1]*np.conj(betas[f2])*DHere #(beta)(beta^*)(D)
        f1Col, f2Col = freq2Col[f1], freq2Col[f2]
        
        rowIndices[termsPerRow*n:termsPerRow*n+termsPerRow/2] = n #the first termsPerRow/2 terms are on this row
        rowIndices[termsPerRow*n+termsPerRow/2:termsPerRow*n+termsPerRow] = n+nPairs #the next termsPerRow/2 terms are on the corresponding imaginary part row

        #Frequency terms
        for i,colIndex in enumerate([f1Col,f2Col,f1Col+1,f2Col+1]): 
            colIndices[termsPerRow*n + i] = colIndex
            colIndices[termsPerRow*n + termsPerRow/2 + i] = colIndex
        for i,coeff in enumerate([np.real(bbstarSD), np.real(bbstarSD), -np.imag(bbstarSD), np.imag(bbstarSD)]): coeffs[termsPerRow*n+i] = coeff
        for i,coeff in enumerate([np.imag(bbstarSD), np.imag(bbstarSD), np.real(bbstarSD), -np.real(bbstarSD)]): coeffs[termsPerRow*n + termsPerRow/2 + i] = coeff        
        
        #Sigma terms
        for i,colIndex in enumerate(realSigmaPolyCols):
            colIndices[termsPerRow*n + 4 + i] = colIndex
            colIndices[termsPerRow*n + termsPerRow/2 + 4 + i] = colIndex
        for i,order in enumerate(range(realSigmaPolyOrder,-1,-1)):
            coeffs[termsPerRow*n + 4 + i] = (u/uMean)**order * np.real(bbstarD)
            coeffs[termsPerRow*n + termsPerRow/2 + 4 + i] = (u/uMean)**order * np.imag(bbstarD)
        for i,colIndex in enumerate(imagSigmaPolyCols):
            colIndices[termsPerRow*n + 4 + realSigmaPolyOrder+1  + i] = colIndex
            colIndices[termsPerRow*n + termsPerRow/2 + 4 + realSigmaPolyOrder+1 + i] = colIndex
        for i,order in enumerate(range(imagSigmaPolyOrder,-1,-1)):
            coeffs[termsPerRow*n + 4 + realSigmaPolyOrder+1  + i] = -(u/uMean)**order * np.imag(bbstarD)
            coeffs[termsPerRow*n + termsPerRow/2 + 4 + realSigmaPolyOrder+1 + i] = (u/uMean)**order * np.real(bbstarD)        
            
        #D terms
        for i,colIndex in enumerate(realDsPolyCols):
            colIndices[termsPerRow*n + 4 + realSigmaPolyOrder + imagSigmaPolyOrder + 2 + i] = colIndex
            colIndices[termsPerRow*n + termsPerRow/2 + 4 + realSigmaPolyOrder + imagSigmaPolyOrder + 2 + i] = colIndex
        for i,order in enumerate(range(realDsPolyOrder,-1,-1)):
            coeffs[termsPerRow*n + 4 + realSigmaPolyOrder + imagSigmaPolyOrder + 2 + i] = (du/duMean)**order * np.real(bbstarS)
            coeffs[termsPerRow*n + termsPerRow/2 + 4 + realSigmaPolyOrder + imagSigmaPolyOrder + 2 + i] = (du/duMean)**order * np.imag(bbstarS)
        for i,colIndex in enumerate(imagDsPolyCols):
            colIndices[termsPerRow*n + 4 + realSigmaPolyOrder + imagSigmaPolyOrder + realDsPolyOrder + 3 + i] = colIndex
            colIndices[termsPerRow*n + termsPerRow/2 + 4 + realSigmaPolyOrder + imagSigmaPolyOrder + realDsPolyOrder + 3 + i] = colIndex
        for i,order in enumerate(range(imagDsPolyOrder,-1,-1)):
            coeffs[termsPerRow*n + 4 + realSigmaPolyOrder + imagSigmaPolyOrder + realDsPolyOrder + 3  + i] = -(du/duMean)**order * np.imag(bbstarS)
            coeffs[termsPerRow*n + termsPerRow/2 + 4 + realSigmaPolyOrder + imagSigmaPolyOrder + realDsPolyOrder + 3 + i] = (du/duMean)**order * np.real(bbstarS) 
                        
    A = csr_matrix((coeffs,(rowIndices,colIndices)))
    Ninv = csr_matrix((np.append(noiseCovDiag**-1,noiseCovDiag**-1), (np.arange(2*nPairs), np.arange(2*nPairs))))
    AtNinvA = (A.conjugate().transpose().dot(Ninv)).dot(A).toarray()    
    if iteration == 0: print 'Zero Eigenvalues: ', len(AtNinvA) - np.linalg.matrix_rank(AtNinvA)
    deltasSplitReIm = np.append(np.real(deltas),np.imag(deltas))
    xHat = np.linalg.pinv(AtNinvA).dot(A.T.conjugate().dot(Ninv.dot(deltasSplitReIm)))
    


    newBetas = {freq: betas[freq]*(1+alpha*(xHat[freq2Col[freq]] + 1.0j*xHat[freq2Col[freq]+1])) for freq in allFreqs}
    betas = newBetas
    realSigmas = realSigmas + alpha*xHat[realSigmaPolyCols]
    imagSigmas = imagSigmas + alpha*xHat[imagSigmaPolyCols]
    realDs = realDs + alpha*xHat[realDsPolyCols]
    imagDs = imagDs + alpha*xHat[imagDsPolyCols]
    
    #Renormalize
    bandpassMean = np.average(np.asarray([np.abs(betas[f]) for f in allFreqs]))    
    for f in allFreqs: betas[f] = betas[f] / bandpassMean
    bandpassAngleMean = np.average(np.asarray([np.angle(betas[f]) for f in allFreqs]))
    for f in allFreqs: betas[f] = betas[f] * np.exp(-1.0j*bandpassAngleMean)
    Dmean = np.average(np.abs(np.asarray([(np.polyval(realDs,entry['du']/duMean) + 1.0j*np.polyval(imagDs,entry['du']/duMean)) for ((f1,f2,u),entry) in baselineFreqPairs.items()])))
    realSigmas *= bandpassMean**2 * Dmean
    imagSigmas *= bandpassMean**2 * Dmean
    realDs /= Dmean
    imagDs /= Dmean

    betaBetaStarList = np.asarray([betas[f1]*np.conj(betas[f2]) for ((f1,f2,u),entry) in baselineFreqPairs.items()])    
    SigmasList = np.asarray([(np.polyval(realSigmas,u/uMean) + 1.0j*np.polyval(imagSigmas,u/uMean)) for ((f1,f2,u),entry) in baselineFreqPairs.items()])
    DsList = np.asarray([(np.polyval(realDs,entry['du']/duMean) + 1.0j*np.polyval(imagDs,entry['du']/duMean)) for ((f1,f2,u),entry) in baselineFreqPairs.items()])
    allChiSq = np.abs(visCorrs - betaBetaStarList * SigmasList * DsList)**2 / noiseCovDiag
    plt.figure(30); plt.clf()
    plt.plot(allFreqs,[np.abs(betas[freq]) for freq in allFreqs],'.')
    plt.plot(allFreqs,[np.abs(binnedLincalBetas[freq]) for freq in allFreqs],'r.')
    plt.legend(['Polynomial Lincal','Discrete Lincal'])
    plt.title('Polynomial Lincal abs(bandpass)')
    plt.xlabel('f (GHz)')
    plt.figure(31); plt.clf()
    plt.plot(np.arange(uMin,uMax,1),np.abs(np.polyval(realSigmas,np.arange(uMin,uMax,1)/uMean) + (np.polyval(imagSigmas,np.arange(uMin,uMax,1)/uMean))))
    plt.title('Polynomial abs(Sigma(u))')
    plt.xlabel('u')
    plt.show() 
    plt.draw()
    plt.pause(0.05)

    print str(iteration) + ') mean chi^2 = ' + str(np.average(allChiSq))
    

#%% Diagnostic plotting

betaBetaStarList = np.asarray([betas[f1]*np.conj(betas[f2]) for ((f1,f2,u),entry) in baselineFreqPairs.items()])    
SigmasList = np.asarray([(np.polyval(realSigmas,u/uMean) + 1.0j*np.polyval(imagSigmas,u/uMean)) for ((f1,f2,u),entry) in baselineFreqPairs.items()])
DsList = np.asarray([(np.polyval(realDs,entry['du']/duMean) + 1.0j*np.polyval(imagDs,entry['du']/duMean)) for ((f1,f2,u),entry) in baselineFreqPairs.items()])
predictedCorrs =  betaBetaStarList * SigmasList * DsList

#use binned version instead:
#predictedCorrs = np.asarray([binnedLincalBetas[f1]*np.conj(binnedLincalBetas[f2])*binnedLincalSigmas[entry['uBin']]*binnedLincalDs[entry['duBin']] for ((f1,f2,u),entry) in baselineFreqPairs.items()]) 
#AtNinvA = binnedLincalAtNinvA


errorList = np.abs(visCorrs - predictedCorrs)
fracErrorList = errorList / np.abs(visCorrs)

f1List = np.asarray([f1 for ((f1,f2,u),entry) in baselineFreqPairs.items()])
f2List = np.asarray([f2 for ((f1,f2,u),entry) in baselineFreqPairs.items()])
uList = np.asarray([u for ((f1,f2,u),entry) in baselineFreqPairs.items()])
uBinList = np.asarray([entry['uBin'] for ((f1,f2,u),entry) in baselineFreqPairs.items()])
duList = np.asarray([entry['du'] for ((f1,f2,u),entry) in baselineFreqPairs.items()])
b1List = np.asarray([entry['separations'][0] for ((f1,f2,u),entry) in baselineFreqPairs.items()])
b2List = np.asarray([entry['separations'][1] for ((f1,f2,u),entry) in baselineFreqPairs.items()])
deltaProdList = np.asarray([entry['averageDeltaProd'] for ((f1,f2,u),entry) in baselineFreqPairs.items()])
noiseList = np.asarray(noiseCovDiag)
#SigmaList = np.asarray([(np.polyval(realSigmas,u/uMean) + 1.0j*np.polyval(imagSigmas,u/uMean)) for ((f1,f2,u),entry) in baselineFreqPairs.items()])

def lincalScatter(x,y,color=None, figNum=100, xs='log', ys='log', title='', clear=True, xl='', yl=''):
    plt.figure(figNum); 
    if clear: plt.clf()
    if color is not None: plt.scatter(x,y,c=color)
    else: plt.scatter(x,y)
    plt.yscale(ys)
    plt.xscale(xs)
    plt.ylim([.9*np.min(y), 1.1*np.max(y)])
    plt.xlim([.9*np.min(x), 1.1*np.max(x)])
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.title(title)

#lincalScatter(uList,allChiSq, figNum=100, color=b1List, xs='linear', title = 'color = b1', xl='u', yl='$\chi^2$')
#lincalScatter(uList,allChiSq, figNum=101, color=b2List, xs='linear', title = 'ChiSq vs. u (color = b2)')
#lincalScatter(duList,np.abs(visCorrs - predictedCorrs), figNum=102, color=uList, xs='linear', title = 'ChiSq vs. du/u (color = u)')
#lincalScatter(noiseList**.5,allChiSq, figNum=103, color=uList, xs='log', title = 'ChiSq vs. Noise Std (color = u)')
#lincalScatter(b1List,allChiSq, figNum=104, color=f2List, xs='linear', title = 'ChiSq vs. b1 (color = f2)')
#lincalScatter(b2List,allChiSq, figNum=105, color=f1List, xs='linear', title = 'ChiSq vs. b2 (color = f1)')
#lincalScatter(SigmaList,allChiSq, figNum=106, color=f1List, xs='log', title = 'ChiSq vs. Sigma (color = f1)')
lincalScatter(deltaProdList,noiseList, figNum=107, color=uList, xs='log', title = 'ChiSq vs. Sigma (color = f1)')


lincalScatter(np.abs(predictedCorrs), np.abs(visCorrs), color=np.log10(np.abs(noiseList)), figNum=108, ys='log', xs='log', title = 'Color is du')
plt.xlabel('predicted visibility correlation')
plt.ylabel('observed visibility correlation')
plt.colorbar()
plt.plot([0,1],[0,1])    

lincalScatter(uList,fracErrorList, figNum=109, color=duList, xs='linear', title = 'Fractional Error vs. u (color = du)')
#TODO: investigate weightings by u or sigma. Also look at what's going on at very low u for high du.



#% Plot Expected Error vs. Observed Std

errorList = np.real(predictedCorrs - visCorrs)
#relativeErrorList = np.abs(predictedCorrs - visCorrs)/np.abs(predictedCorrs)
freqCompiledList = {freq: [] for freq in allFreqs}
for f1,f2,error in zip(f1List,f2List,errorList):
    freqCompiledList[f1].append(error)
    freqCompiledList[f2].append(error)    

freqAvgErrors = np.asarray([np.mean(np.asarray(freqCompiledList[freq])) for freq in allFreqs])
freqMedianErrors = np.asarray([np.median(np.asarray(freqCompiledList[freq])) for freq in allFreqs])
freqStdErrors = np.asarray([np.std(np.asarray(freqCompiledList[freq])) for freq in allFreqs])
freqCounts = [len(freqCompiledList[freq]) for freq in allFreqs]

plt.figure(1123); plt.clf()
#plt.errorbar(allFreqs,freqAvgErrors,yerr=freqStdErrors)
plt.plot(allFreqs,freqAvgErrors,'k.')
plt.plot(allFreqs,freqMedianErrors,'g.-')
plt.plot(allFreqs,freqStdErrors,'b-')

plt.plot(allFreqs,.03*np.asarray([np.abs(betas[freq]) for freq in allFreqs])*((np.diag(np.linalg.pinv(AtNinvA))[0:2*nFreqs:2]))**.5,'r-')
plt.legend(['Average Real Error','Median Real Error','Std of Real Error','Rescaled Expected Noise'])
plt.xlabel('f (GHz)')

#%% Examine error correlations

plt.figure(99); plt.clf()
AtNinvAinv = np.linalg.pinv(AtNinvA)[0:2*nFreqs:2,0:2*nFreqs:2]
inverseSqrtDiag = np.diag(np.diag(AtNinvAinv)**-.5)
plt.imshow(inverseSqrtDiag.dot(AtNinvAinv.dot(inverseSqrtDiag)), interpolation='none', extent = [allFreqs[0],allFreqs[-1],allFreqs[0],allFreqs[-1]])
plt.title('Error Correlation Matrix for Real Part')
plt.colorbar()


#%% Find worst frequency channels

freq2chan = {freq: [] for freq in allFreqs}
for ((f1,f2,u),entry) in baselineFreqPairs.items():    
    freq2chan[f1]=(entry['chans'][0])
    freq2chan[f2]=(entry['chans'][1])
#print [freq2chan[freq] for freq in allFreqs[freqAvgErrors>7e-5]]

##%% Testing Real/Imag vs. Abs/Angle
#
#plt.figure(1212); plt.clf()
##plt.plot(duBinCenters.values(),np.abs(np.asarray([binnedLincalDs[duBin] for duBin in allduBins])),'r.')
#plt.plot(duBinCenters.values(),np.real(np.asarray([binnedLincalDs[duBin] for duBin in allduBins])),'b-')
#plt.plot(duBinCenters.values(),np.angle(np.asarray([binnedLincalDs[duBin] for duBin in allduBins])),'b-')
#plt.plot(duBinCenters.values(),np.imag(np.asarray([binnedLincalDs[duBin] for duBin in allduBins])),'g-')
#


#%%
deltaProdFreqAvg = {freq: [] for freq in allFreqs}
for ((f1,f2,u),entry) in baselineFreqPairs.items(): 
    if entry['du'] < .5:
        deltaProdFreqAvg[f1].append(entry['averageDeltaProd'])
        deltaProdFreqAvg[f2].append(entry['averageDeltaProd'])
deltaProdFreqAvg = [np.average(np.asarray(deltaProdFreqAvg[freq])) for freq in allFreqs]
plt.figure(332); plt.clf()
plt.plot(allFreqs,deltaProdFreqAvg,'.')
plt.plot(allFreqs,.05*np.asarray([np.abs(betas[freq]) for freq in allFreqs])*((np.diag(np.linalg.pinv(AtNinvA))[0:2*nFreqs:2]))**.5,'r-')