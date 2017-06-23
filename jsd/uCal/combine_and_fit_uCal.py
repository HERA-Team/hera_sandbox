#! /usr/bin/env python
import numpy as np
import capo.uCal as uc
import matplotlib.pyplot as plt
import optparse, sys, os

o = optparse.OptionParser()
o.set_usage('combine_and_fit_uCal.py [.uCalResults.npz files] [options]')
o.add_option('--nofit', action='store_true',
    help='For speed, does not try to fit the bandpass.')
o.add_option('--combinedSuffix', type='string', default='combined',
    help='Saves the result as the first uCalResults.npz file, but with .npz replaced by .[suffix].npz. Default "combined"')
opts,args = o.parse_args(sys.argv[1:])

if len(args) == 0:
    print "WARNING!!! No arguments provided. Using testing parameters."
    resultsFiles = ["./Data/" + file for file in os.listdir("./Data") if file.startswith('zen.2456679.3') and file.endswith('xx.uCalResults.npz')]
else: resultsFiles = args

#%% Combine per-integration results
allSourceFiles, allBandpasses, meanBandpass, stdBandpass, unflaggedChans, channelRMSs, overallChannelRMS, bandpassFit = uc.from_npz(resultsFiles[0])
for file in resultsFiles[1:]:
    sourceFiles, allBP, _, _, unflagged, RMSs, _, _ = uc.from_npz(file)
    for source in sourceFiles: allSourceFiles = np.append(allSourceFiles,source)
    allBandpasses = np.concatenate((allBandpasses,allBP))
    unflaggedChans = np.concatenate((unflaggedChans,unflagged))
    channelRMSs = np.concatenate((channelRMSs,RMSs))

#%% Recompute Summary Results
meanBandpass = np.mean(np.ma.masked_equal(allBandpasses,0),axis=0).data
stdBandpass = np.std(np.abs(np.ma.masked_invalid(np.ma.masked_equal(allBandpasses,0))),axis=0).data
overallChannelRMS = np.mean(np.ma.masked_invalid(np.ma.masked_equal(channelRMSs,0))**2,axis=0).data ** .5

#%% Perform Fit
bandpassFit = None
if not opts.nofit:
    print "Fitting not yet enabled!!!"
    #TODO: enable fitting

#%% Save Results
outfile = resultsFiles[0].replace('.npz', '.' + opts.combinedSuffix + '.npz')
print 'Combined the results of ' + str(len(allSourceFiles)) + ' files with ' + str(allBandpasses.shape[0]) + ' integrations.'
uc.to_npz(outfile, allSourceFiles, allBandpasses, meanBandpass, stdBandpass, unflaggedChans, channelRMSs, overallChannelRMS, bandpassFit)




# import cPickle as pickle
# import matplotlib.pyplot as plt
# import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
# import scipy
# import scipy.constants
# import aipy as a
# import optparse, sys, os
# import capo
# import capo.uCal as uc
# from joblib import Parallel, delayed
# import multiprocessing
# from scipy.optimize import curve_fit

# try: 
#     resultsFile = sys.argv[1]
#     if resultsFile.endswith('.uCalResults.npz'): resultsFile = resultsFile.replace('.uCalResults.npz','diagnosticResults.p')
#     elif resultsFile.endswith('.npz'): resultsFile = resultsFile.replace('.npz','diagnosticResults.p')
#     results = pickle.load(open(resultsFile,'rb'))
# except: 
#     results = pickle.load(open('./Data/zen.2456943.38268.xx.uvcRRE.diagnosticResults.p','rb'))
# #    results = pickle.load(open('./Data/zen.2456679.35658.xx.diagnosticResults.p','rb'))
# for key,val in results.items(): exec(key + '=val')

# #%% Plot fits
# plt.figure(755); plt.clf()
# plt.plot(degreeRange,BICs)
# plt.xlabel('Degree'); plt.ylabel('BIC')
# plt.title('Bayesian Information Criterion for Different Numbers of Modes Fit')

# plt.figure(766); plt.clf()
# plt.errorbar(uCal.chans, np.real(betas)*freqs[uCal.chans]**2.55, yerr=observedRealErrors*freqs[uCal.chans]**2.55)
# plt.errorbar(uCal.chans, np.imag(betas*freqs[uCal.chans]**2.55),yerr=observedImagErrors*freqs[uCal.chans]**2.55)
# plt.plot(uCal.chans, np.real(model),'r')
# plt.plot(uCal.chans, np.imag(model),'k')
# plt.title('Cosine Fit to Bandpass')
# plt.xlabel('channel'); plt.ylabel('beta')
# plt.legend(['Re(beta) fit', 'Im(beta) fit','Re(beta)','Im(beta)'])

# def lincalScatter(x,y,color=None, figNum=100, xs='log', ys='log', title='', clear=True, xl='', yl=''):
#     plt.figure(figNum); 
#     if clear: plt.clf()
#     if color is not None: plt.scatter(x,y,c=color)
#     else: plt.scatter(x,y)
#     plt.yscale(ys); plt.xscale(xs)
#     plt.ylim([.9*np.min(y), 1.1*np.max(y)]); plt.xlim([.9*np.min(x), 1.1*np.max(x)])
#     plt.xlabel(xl); plt.ylabel(yl); plt.title(title)
# duList = np.asarray([entry['du'] for entry in uCal.blChanPairs.values()])
# uList = np.asarray([entry['u'] for entry in uCal.blChanPairs.values()])
# ch1List = np.asarray([ch1 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])
# ch2List = np.asarray([ch2 for (ch1,bl1,ch2,bl2) in uCal.blChanPairs.keys()])

# #%%Betas
# plt.figure(102); plt.clf()
# #inferredErrorsOnAbsBeta = np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
# #plt.errorbar(uCal.chans,np.abs(betas),yerr=inferredErrorsOnAbsBeta)
# plt.plot(uCal.chans,np.angle(betas))
# plt.plot(uCal.chans, np.angle(np.asarray(bootstrapBetas).T),':.')#'--.')
# plt.xlabel('chan')
# plt.ylabel('Angle(Beta)'); plt.title('Beta with Bootstraps')

# #%%Sigmas
# plt.figure(103); plt.clf()
# inferredErrorsOnAbsSigma = ((np.diag(np.linalg.pinv(uCal.AtNinvA))[2*uCal.nChans : 2*uCal.nChans+2*uCal.nuBins : 2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1+2*uCal.nChans : 2*uCal.nChans+2*uCal.nuBins : 2]))**.5
# plt.errorbar(np.asarray(uCal.uBins)[:,0], np.abs(Sigmas),yerr=inferredErrorsOnAbsSigma)
# plt.plot(np.asarray(uCal.uBins)[:,0], np.abs(np.asarray(bootstrapSigmas).T),':.')#'--.')
# plt.xlabel('uBin')
# plt.ylabel('Abs(Sigma)'); plt.title('Sigma with Bootstraps')     

# #%%Ds
# plt.figure(104); plt.clf()
# inferredErrorsOnAbsD = ((np.diag(np.linalg.pinv(uCal.AtNinvA))[2*uCal.nChans+2*uCal.nuBins :: 2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1+2*uCal.nChans+2*uCal.nuBins :: 2]))**.5
# plt.plot(uCal.duBins, np.abs(np.asarray(bootstrapDs).T),':.')#'--.')
# plt.xlabel('duBin')
# plt.ylabel('Abs(Ds)'); plt.title('Ds with Bootstraps')   

# #%%Bandpass
# plt.figure(1); plt.clf()
# inferredErrorsOnBandpass = freqs[uCal.chans]**2.55 * np.abs(betas)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
# #plt.errorbar(uCal.chans,np.abs(freqs[uCal.chans]**2.55 * betas) / np.mean(np.abs(freqs[uCal.chans]**2.55 * betas)),yerr=inferredErrorsOnBandpass/np.mean(np.abs(freqs[uCal.chans]**2.55 * betas)))
# plt.errorbar(uCal.chans,np.abs(model) / np.mean(np.abs(model)),yerr=inferredErrorsOnBandpass/np.mean(np.abs(model)))
# plt.xlabel('Channel'); 
# plt.ylabel('Abs(Lincal Bandpass)'); plt.title(r'Lincal Bandpass Corrected by $\nu^{2.55}$')

# #%%Predicted vs. Observed Scatter
# plt.figure(2); plt.clf()
# visCorrs = np.asarray([entry['visCorr'] for entry in uCal.blChanPairs.values()])
# predictedCorrs = visCorrs - uCal.computeErrors(betas, Sigmas, Ds)
# lincalScatter(np.abs(predictedCorrs), np.abs(visCorrs), color=np.linalg.norm(uList,axis=1), figNum=2, ys='log', xs='log', title = 'Scatter in Observed/Modeled Visibility Correlations')
# plt.plot([0,1],[0,1],'k--')
# plt.xlabel('Abs(Predicted Correlations)'); plt.ylabel('Abs(Observe Correlations)')
# plt.colorbar(label='u (wavelengths)')
# #

# #%%Examine error correlations
# plt.figure(3); plt.clf()
# AtNinvAinv = np.linalg.pinv(uCal.AtNinvA)[0:2*uCal.nChans:2,0:2*uCal.nChans:2]
# inverseSqrtDiag = np.diag(np.diag(AtNinvAinv)**-.5)
# plt.imshow(inverseSqrtDiag.dot(AtNinvAinv.dot(inverseSqrtDiag)), interpolation='nearest', extent = [uCal.chans[0],uCal.chans[-1],uCal.chans[0],uCal.chans[-1]], vmin=-.2, vmax=1)
# plt.title('Error Correlation Matrix for Real Part of Beta')
# plt.xlabel('Channel'); plt.ylabel('Channel')
# plt.colorbar()

# #%%Channel averaged Chi^2
# plt.figure(1000); plt.clf()
# plt.plot(uCal.chans, uCal.chanAvgErrors,'.-')
# plt.title('Average Chi^2 By Channel')
# plt.xlabel('Channel'); plt.ylabel('Avg. Chi^2')
# plt.show()


# #%% Compare different paremters effects on the Bandpass
# if False:
#     plt.figure(222); plt.clf()
#     #prefixes = ['umin0', 'umin10', 'umin20', 'umin30']
#     #prefixes = ['umax100', 'umax120', 'umax140', 'umax160']
#     prefixes = ['deltaumax.2', 'deltaumax.3', 'deltaumax.4', 'deltaumax.5', 'deltaumax.6', 'deltaumax.7']
#     for prefix in prefixes:
#         results = pickle.load(open('./Data/zen.2456679.35658.xx.'+prefix+'.diagnosticResults.p','rb'))
#         for key,val in results.items(): exec(key + '=val')
#         betasRescaled = betas / np.mean(np.abs(betas))
#         inferredErrorsOnBandpass = freqs[uCal.chans]**2.55 * np.abs(betasRescaled)*((np.diag(np.linalg.pinv(uCal.AtNinvA))[0:2*uCal.nChans:2]) + (np.diag(np.linalg.pinv(uCal.AtNinvA))[1:2*uCal.nChans:2]))**.5
#         plt.errorbar(uCal.chans,np.abs(freqs[uCal.chans]**2.55 * betasRescaled),yerr=inferredErrorsOnBandpass)
#     plt.xlabel('Channel'); 
#     plt.ylabel('Abs(Lincal Bandpass)'); plt.title(r'Lincal Bandpass Corrected by $\nu^{2.55}$')
#     plt.legend(prefixes)