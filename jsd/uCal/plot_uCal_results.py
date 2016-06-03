#%% Plot fits
plt.figure(755); plt.clf()
plt.plot(degreeRange,BICs)
plt.xlabel('Degree'); plt.ylabel('BIC')

plt.figure(766); plt.clf()
plt.errorbar(uCal.chans, np.real(betas)*freqs[uCal.chans]**2.55, yerr=observedRealErrors*freqs[uCal.chans]**2.55)
plt.errorbar(uCal.chans, np.imag(betas*freqs[uCal.chans]**2.55),yerr=observedImagErrors*freqs[uCal.chans]**2.55)
plt.plot(uCal.chans, np.real(model),'r')
plt.plot(uCal.chans, np.imag(model),'k')
plt.xlabel('channel'); plt.ylabel('beta')
plt.legend(['Re(beta) fit', 'Im(beta) fit','Re(beta)','Im(beta)'])
#
#%%##########################################
#   Diagnostic Plotting
#############################################
if True:

    def lincalScatter(x,y,color=None, figNum=100, xs='log', ys='log', title='', clear=True, xl='', yl=''):
        plt.figure(figNum); 
        if clear: plt.clf()
        if color is not None: plt.scatter(x,y,c=color)
        else: plt.scatter(x,y)
        plt.yscale(ys); plt.xscale(xs)
        plt.ylim([.9*np.min(y), 1.1*np.max(y)]); plt.xlim([.9*np.min(x), 1.1*np.max(x)])
        plt.xlabel(xl); plt.ylabel(yl); plt.title(title)
    duList = np.asarray([entry['du'] for entry in uCal.blChanPairs.values()])
    uList = np.asarray([entry['u'] for entry in uCal.blChanPairs.values()])
    ch1List = np.asarray([ch1 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
    ch2List = np.asarray([ch2 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
#     relativeScatter = np.zeros(len(freqs))
#     relativeScatter[uCal.chans] = observedVarianceByChannel/modelVarianceByChannel
#     relativeScatterList = np.asarray([relativeScatter[ch1] + relativeScatter[ch2] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])

    #%%Betas
    plt.figure(102); plt.clf()
    inferredErrorsOnAbsBeta = np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
    plt.errorbar(uCal.chans,np.abs(betas),yerr=inferredErrorsOnAbsBeta)
    plt.plot(uCal.chans, np.abs(np.asarray(bootstrapBetas).T),':.')#'--.')
    plt.xlabel('chan')#     plt.xlabel('Frequency (GHz)'); 
    plt.ylabel('Abs(Beta)'); plt.title('Beta with Bootstraps')

   #%%Sigmas
    plt.figure(103); plt.clf()
    inferredErrorsOnAbsSigma = ((np.diag(np.linalg.pinv(uCal.AtNinvA))[2*uCal.nChans : 2*uCal.nChans+2*uCal.nuBins : 2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1+2*uCal.nChans : 2*uCal.nChans+2*uCal.nuBins : 2]))**.5
    plt.errorbar(np.asarray(uCal.uBins)[:,0], np.abs(Sigmas),yerr=inferredErrorsOnAbsSigma)
    plt.plot(np.asarray(uCal.uBins)[:,0], np.abs(np.asarray(bootstrapSigmas).T),':.')#'--.')
    plt.xlabel('uBin')#     plt.xlabel('Frequency (GHz)'); 
    plt.ylabel('Abs(Sigma)'); plt.title('Sigma with Bootstraps')     

   #%%Ds
    plt.figure(104); plt.clf()
    inferredErrorsOnAbsD = ((np.diag(np.linalg.pinv(uCal.AtNinvA))[2*uCal.nChans+2*uCal.nuBins :: 2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1+2*uCal.nChans+2*uCal.nuBins :: 2]))**.5
#    plt.errorbar(uCal.duBins, np.abs(Ds),yerr=inferredErrorsOnAbsD)
    plt.plot(uCal.duBins, np.abs(np.asarray(bootstrapDs).T),':.')#'--.')
#    plt.errorbar(uCal.duBins, np.abs(allResults[-1][2]), yerr=inferredErrorsOnAbsD)
#    plt.plot(uCal.duBins, np.abs(np.asarray(allResults[-1][6]).T),':.')#'--.')
    plt.xlabel('uBin')#     plt.xlabel('Frequency (GHz)'); 
    plt.ylabel('Abs(Ds)'); plt.title('Ds with Bootstraps')   

    #%%Bandpass
    plt.figure(1); plt.clf()
    inferredErrorsOnBandpass = freqs[uCal.chans]**2.55 * np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
    plt.errorbar(np.arange(.1,.2,.1/203)[uCal.chans],np.abs(freqs[uCal.chans]**2.55 * betas),yerr=inferredErrorsOnBandpass)
#     plt.errorbar(uCal.chans,np.abs(freqs[uCal.chans]**2.55 * betas),yerr=inferredErrorsOnAbsBeta)
#     plt.plot(np.arange(.1,.2,.1/203)[uCal.chans], np.asarray([freqs[uCal.chans]**2.55 * beta for beta in np.abs(np.asarray(splitBetas))]).T,'.')
    plt.xlabel('Frequency (GHz)'); 
    plt.ylabel('Abs(Lincal Bandpass)'); plt.title(r'Lincal Bandpass Corrected by $\nu^{2.55}$')

    #%%Predicted vs. Observed Scatter
    plt.figure(2); plt.clf()
    visCorrs = np.asarray([entry['visCorr'] for entry in uCal.blChanPairs.values()])
    predictedCorrs = visCorrs - uCal.computeErrors(betas, Sigmas, Ds)
    lincalScatter(np.abs(predictedCorrs), np.abs(visCorrs), color=np.linalg.norm(uList,axis=1), figNum=2, ys='log', xs='log', title = 'color = u') #color=np.log10(relativeScatterList)
    plt.plot([0,1],[0,1],'k--')
    plt.xlabel('Abs(Predicted Correlations)'); plt.ylabel('Abs(Observe Correlations)')
#

    #%%Examine error correlations
    plt.figure(3); plt.clf()
    AtNinvAinv = np.linalg.pinv(uCal.AtNinvA)[0:2*uCal.nChans:2,0:2*uCal.nChans:2]
    inverseSqrtDiag = np.diag(np.diag(AtNinvAinv)**-.5)
    plt.imshow(inverseSqrtDiag.dot(AtNinvAinv.dot(inverseSqrtDiag)), interpolation='nearest', extent = [uCal.chans[0],uCal.chans[-1],uCal.chans[0],uCal.chans[-1]], vmin=-.2, vmax=1)
    plt.title('Error Correlation Matrix for Real Part of Beta')
    plt.xlabel('Channel'); plt.ylabel('Channel')
    plt.colorbar()
    
    #%%
#    plt.figure(4); plt.clf()
#    plt.plot(uCal.chans, inferredErrorsOnAbsBeta)
#    plt.plot(uCal.chans, np.std(np.abs(np.asarray(bootstrapBetas).T),axis=1),'.-')
#    plt.legend(['Inferred Errors','Bootstrap Stds'])
#    plt.xlabel('Channel')
#
#    plt.figure(5); plt.clf()
#    plt.plot(np.asarray(uCal.uBins)[:,0], inferredErrorsOnAbsSigma)
#    plt.plot(np.asarray(uCal.uBins)[:,0], np.std(np.abs(np.asarray(bootstrapSigmas).T),axis=1),'.-')
#    plt.legend(['Inferred Errors','Bootstrap Stds'])
#    plt.xlabel('uBin')
#
#    plt.figure(6); plt.clf()
#    plt.plot(uCal.duBins, inferredErrorsOnAbsD)
#    plt.plot(uCal.duBins, np.std(np.abs(np.asarray(bootstrapDs).T),axis=1),'.-')
#    plt.legend(['Inferred Errors','Bootstrap Stds'])
#    plt.xlabel('duBin')
#    
#    
#    #%% Compare errors to various other quantities
#    
#    #%% Error vs. Noise
#    plt.figure(25); plt.clf()
#    lincalScatter(noiseCovDiag**.5, np.abs(uCal.computeErrors(betas, Sigmas, Ds)), figNum=25, color=np.linalg.norm(uList,axis=1))
#    plt.plot([0,1],[0,1],'k--')
#    plt.xlabel('Expected Std'); plt.ylabel('Abs(Error on VisCorr)'); plt.title('Noise Model Based on Bootstrapped Variances (color = u)')
#    
#    plt.figure(26); plt.clf()
#    lincalScatter(uCal.generateNoiseCovariance(betas)**.5, np.abs(uCal.computeErrors(betas, Sigmas, Ds)), figNum=26, color=np.linalg.norm(uList,axis=1))
#    plt.plot([0,1],[0,1],'k--')
#    plt.xlabel('Expected Std'); plt.ylabel('Abs(Error on VisCorr)'); plt.title('Noise Model Based on Beta and T_obs (color = u)')
#
#    #%% Error vs. Data
#    plt.figure(27); plt.clf()
#    lincalScatter(np.abs(visCorrs), np.abs(uCal.computeErrors(betas, Sigmas, Ds)), color=np.linalg.norm(uList,axis=1), figNum=27)
#    plt.xlabel('abs(visibility correlation)'); plt.ylabel('abs(error)'); plt.title('Bootstrapped Variance Model (color = u)')  
#    plt.plot([0,1],[0,1],'k--')    
#    
#    plt.figure(28); plt.clf()
#    lincalScatter(np.abs(visCorrs), np.abs(uCal.computeErrors(allResults[0][0], allResults[0][1], allResults[0][2])), color=np.linalg.norm(uList,axis=1), figNum=28)
#    plt.xlabel('abs(visibility correlation)'); plt.ylabel('abs(error)'); plt.title('Noise Model Based on Beta and T_obs (color = u)')
#    plt.plot([0,1],[0,1],'k--')
#    
#    #%% Error vs. Params
#    betaIndex = {uCal.chans[i]: i for i in range(uCal.nChans)}
#    SigmaIndex = {uCal.uBins[i]: i for i in range(uCal.nuBins)}
#    DIndex = {uCal.duBins[i]: i for i in range(uCal.nduBins)}
#
#    trial = 0
#    figNum = 299
#    
#    beta1List = [allResults[trial][0][betaIndex[ch1]] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
#    beta2List = [allResults[trial][0][betaIndex[ch2]] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
#    betaProduct = [allResults[trial][0][betaIndex[ch1]]*np.conj(allResults[trial][0][betaIndex[ch2]]) for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
#    SigmaList = [allResults[trial][1][SigmaIndex[entry['uBin']]] for entry in uCal.blChanPairs.values()]
#    DList = [allResults[trial][2][DIndex[entry['duBin']]] for entry in uCal.blChanPairs.values()]
#    samplesList = [entry['samples'] for entry in uCal.blChanPairs.values()]
#    duList = np.asarray([entry['du'] for entry in uCal.blChanPairs.values()])
#    uList = np.asarray([entry['u'] for entry in uCal.blChanPairs.values()])
#    ch1List = np.asarray([ch1 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
#    ch2List = np.asarray([ch2 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
#    visCorrs = [entry['visCorr'] for entry in uCal.blChanPairs.values()]
#    
#    betaStds = np.std(np.abs(np.asarray(allResults[trial][4]).T),axis=1)
#    SigmaStds = np.std(np.abs(np.asarray(allResults[trial][5]).T),axis=1)
#    DStds = np.std(np.abs(np.asarray(allResults[trial][6]).T),axis=1)
#    beta1StdList = [betaStds[betaIndex[ch1]] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
#    beta2StdList = [betaStds[betaIndex[ch2]] for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()]
#    SigmaStdList = [SigmaStds[SigmaIndex[entry['uBin']]] for entry in uCal.blChanPairs.values()]
#    DStdList = [DStds[DIndex[entry['duBin']]] for entry in uCal.blChanPairs.values()]
#    absErrorList = np.abs(uCal.computeErrors(allResults[trial][0], allResults[trial][1], allResults[trial][2]))
#
#
#
#    def errorComparisonPlot(x,y,color,figNum,xlabel,ylabel,title):
#        plt.figure(figNum); plt.clf()
#        lincalScatter(np.abs(x), np.abs(y), color=color, figNum=figNum)
#        plt.xlabel(xlabel); plt.ylabel(ylabel); plt.title('iteration ' + str(trial) + ': ' + title)
#        plt.plot([0,1e6*np.median(x)],[0,1e6*np.median(y)],'k--')
#        plt.colorbar()
#
##    figNum+=1; errorComparisonPlot(np.asarray(samplesList)**-.5, absErrorList, ch1List, figNum, '1/sqrt(samples)', 'abs(error)', 'color = chan1')
#
##    figNum+=1; errorComparisonPlot(beta1List, absErrorList, ch1List, figNum, 'beta1', 'abs(error)', 'color = chan1')
##    figNum+=1; errorComparisonPlot(beta2List, absErrorList, ch2List, figNum, 'beta2', 'abs(error)', 'color = chan2')
##    figNum+=1; errorComparisonPlot(betaProduct, absErrorList, np.abs(beta1List), figNum, 'beta1*beta2', 'abs(error)', 'color = beta1')
#    figNum+=1; errorComparisonPlot(np.asarray(betaProduct) * np.asarray(samplesList)**-.5, absErrorList, np.abs(beta1List), figNum, 'beta1*beta2/sqrt(samples)', 'abs(error)', 'color = beta1')
#    figNum+=1; errorComparisonPlot(SigmaList, absErrorList, np.linalg.norm(uList,axis=1), figNum, 'Sigma', 'abs(error)', 'color = u')
##    figNum+=1; errorComparisonPlot(DList, absErrorList, duList, figNum, 'Ds', 'abs(error)', 'color = du')
#
##    figNum+=1; errorComparisonPlot(ch1List, absErrorList, beta1List, figNum, 'chan1', 'abs(error)', 'color = beta1')
##    figNum+=1; errorComparisonPlot(ch2List, absErrorList, beta2List, figNum, 'chan2', 'abs(error)', 'color = beta2')
##    figNum+=1; errorComparisonPlot(np.linalg.norm(uList,axis=1), absErrorList, SigmaList, figNum, 'u', 'abs(error)', 'color = Sigma')
##    figNum+=1; errorComparisonPlot(duList, absErrorList, DList, figNum, 'du', 'abs(error)', 'color = D')
#    
##    figNum+=1; errorComparisonPlot(beta1StdList, absErrorList, ch1List, figNum, 'beta1Std', 'abs(error)', 'color = chan1')
##    figNum+=1; errorComparisonPlot(beta1StdList, absErrorList, ch2List, figNum, 'beta1Std', 'abs(error)', 'color = chan2')
#    figNum+=1; errorComparisonPlot(SigmaStdList, absErrorList, np.linalg.norm(uList,axis=1), figNum, 'SigmaStd', 'abs(error)', 'color = u')
##    figNum+=1; errorComparisonPlot(DStdList, absErrorList, duList, figNum, 'DStd', 'abs(error)', 'color = du')
#    
##    figNum+=1; errorComparisonPlot(np.abs(np.asarray(ch1List)-np.asarray(ch2List)), absErrorList, np.linalg.norm(uList,axis=1), figNum, '|delta channel|', 'abs(error)', 'color = u')
#    figNum+=1; errorComparisonPlot(np.abs(visCorrs), np.asarray(absErrorList)/np.asarray(np.abs(visCorrs)), np.linalg.norm(uList,axis=1), figNum, '|visCorr|', '|(error/visCorr)|', 'color = u')
#    figNum+=1; errorComparisonPlot(np.linalg.norm(uList,axis=1), np.asarray(absErrorList)/np.asarray(np.abs(visCorrs)), betaProduct, figNum, '|u|', '|(error/visCorr)|', 'color = betaProduct')
#    figNum+=1; errorComparisonPlot(np.linalg.norm(uList,axis=1), np.asarray(betaProduct) * np.asarray(samplesList)**-.5, np.asarray(absErrorList)/np.asarray(np.abs(visCorrs)), figNum, '|u|', 'beta1*beta2/sqrt(samples)', 'color = |(error/visCorr)|')
#
#    remainingError = np.asarray(absErrorList) / np.abs(np.asarray(betaProduct) * np.asarray(samplesList)**-.5)
#    figNum+=1; errorComparisonPlot(SigmaList, remainingError, np.linalg.norm(uList,axis=1), figNum, 'SigmaList', 'abs(error)/(beta1*beta2*samples**-.5)', 'color = u')
#    figNum+=1; errorComparisonPlot(SigmaStdList, remainingError, np.linalg.norm(uList,axis=1), figNum, 'SigmaStd', 'abs(error)/(beta1*beta2*samples**-.5)', 'color = u')
#    
#
#
#    #%%
#    from scipy.stats import linregress
#    toTest = {'1/sqrt(samples)': np.asarray(samplesList)**-.5,
#            'beta1': np.abs(beta1List),
#            'beta2': np.abs(beta2List),
#            'betaProduct': np.abs(betaProduct),
#            'betaProduct/samples**.5:': np.abs(np.asarray(betaProduct) * np.asarray(samplesList)**-.5),
#            'deltaChan': np.abs(np.asarray(ch1List)-np.asarray(ch2List)),
#            'Sigma': np.abs(SigmaList),
#            'Ds': np.abs(DList),     
#            'beta1Std': beta1StdList,
#            'beta2Std': beta2StdList,
#            'SigmaStd': SigmaStdList,
#            'DStd': DStdList,
#            'beta1Std*beta2Std': np.asarray(beta1StdList)*np.asarray(beta2StdList),
#            'betaProduct*abs(Sigma)/samples**.5':  np.abs(np.asarray(betaProduct) * np.asarray(samplesList)**-.5) / np.asarray(np.abs(SigmaList)),
#            'abs(visCorr)': np.abs(visCorrs),
#            'u': np.linalg.norm(uList,axis=1),
#            }
#
#    #remainingError = absErrorList
##    remainingError = np.asarray(absErrorList) / np.abs(np.asarray(betaProduct) * np.asarray(samplesList)**-.5)
#    remainingError = np.asarray(absErrorList)/np.asarray(np.abs(visCorrs))
#
#    for testName, testVector in toTest.items():
#        print '\n', testName
#        slope, intercept, r_value, p_value, std_err = linregress(testVector,remainingError)
#        print 'r =', r_value
#        print 'slope =', slope
#            
#    
#    
#    
#    
#    
##%%
##     inferredBetas = {chan: [] for chan in uCal.chans}
##     betaDict = {uCal.chans[i]: betas[i] for i in range(uCal.nChans)}
##     for (ch1,bl1,ch2,bl2),entry in uCal.blChanPairs.items():
##         inferredBetas[ch1].append(entry['visCorr'] / betaDict[ch2].conj() / SigmaDict[entry['uBin']] / DDict[entry['duBin']])
##         inferredBetas[ch2].append(entry['visCorr'] / betaDict[ch1] / SigmaDict[entry['uBin']] / DDict[entry['duBin']])
##     inferredBetasStd = np.asarray([np.std(inferredBetas[chan])/(len(inferredBetas[chan])**.5) for chan in uCal.chans])
##        plt.plot(uCal.chans, inferredBetasStd)
#
##%%
##chanCompiledList = {chan: [] for chan in uCal.chans}
##chanCompiledList2 = {chan: [] for chan in uCal.chans}
##errorList = uCal.computeErrors(betas, Sigmas, Ds)
##visCorrs = np.asarray([entry['visCorr'] for entry in uCal.blChanPairs.values()])
##modelList = visCorrs# - errorList
##ch1List = np.asarray([ch1 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
##ch2List = np.asarray([ch2 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
##    
##for f1,f2,error,Nii,model in zip(ch1List,ch2List,errorList,noiseCovDiag,modelList):
##    chanCompiledList[f1].append(error/model)#((2*Nii)**.5))
##    chanCompiledList[f2].append(error/model)#((2*Nii)**.5))
##    chanCompiledList2[f1].append(error/(2*Nii)**.5)
##    chanCompiledList2[f2].append(error/(2*Nii)**.5)
##chanAvgErrors = np.asarray([np.mean(np.abs(np.asarray(chanCompiledList[chan]))) for chan in uCal.chans])
##chanAvgErrors2 = np.asarray([np.mean(np.abs(np.asarray(chanCompiledList2[chan])))for chan in uCal.chans])
##plt.figure(4); plt.clf()
##plt.plot(uCal.chans, chanAvgErrors,'.')
##plt.plot(uCal.chans, chanAvgErrors2,'.-')
##plt.ylabel('Channel-Averaged, Noise Weighted Errors'); plt.xlabel('Channel')
###     badChans = np.asarray(uCal.chans)[chanAvgErrors > 2.5]
###     if len(badChans) > 0: print 'Channels with average sigma > 2.5: ', badChans
##
##print np.asarray(uCal.chans)[chanAvgErrors > .5]
##
##%%
##chanCompiledList = {chan: [] for chan in uCal.chans}
##for (ch1,bl1,ch2,bl2),entry in uCal.blChanPairs.items():
##    chanCompiledList[ch1].append(entry['samples'])
##    chanCompiledList[ch2].append(entry['samples'])
##chanTotalSamples = np.asarray([np.sum(np.asarray(chanCompiledList[chan])) for chan in uCal.chans])
##plt.figure(444); plt.clf()
##plt.plot(uCal.chans, chanTotalSamples,'.')
##print np.asarray(uCal.chans)[chanTotalSamples < np.mean(chanTotalSamples)/2]