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