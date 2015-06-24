#! /usr/bin/env python
"""This is a really silly script that's used to make
the necessary directories for the analysis results""" 

import sys,os,numpy as np

data_loc = "/global/homes/a/acliu/globalSig/fq_120_150_testCase"
beam_sigs = [1.57] #(np.pi/18,np.pi/6,5*np.pi/18,7*np.pi/18)
sqGridSideLens = [12] #(4,8,12,16)
variableBeams = [0] #(0,1)
lowerFreq = 120.
upperFreq = 150. #150.
freqSpace = 1.
fqs = np.arange(lowerFreq,upperFreq+freqSpace,freqSpace)
fqs /= 1000. # Convert from MHz to GHz
lowerFreq /= 1000.
mode="diagonal"


if 'plots'.format(data_loc) not in os.listdir(data_loc):
    os.mkdir('{0}/plots'.format(data_loc))

for beam_sig in beam_sigs:
    del_bl = 1/(2*np.pi*beam_sig*lowerFreq)
    for sqGridSideLen in sqGridSideLens:
        for variableBeam in variableBeams:
            #if variableBeam == 0:
            #    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_fixedWidth_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
            #elif variableBeam == 1:
            #    savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_lambdaBeam_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)
            savekey = 'grid_del_bl_{0:.2f}_sqGridSideLen_{1}_beam_sig_{2:.2f}'.format(del_bl,sqGridSideLen,beam_sig)

            # Create directories for plots
            if '{1}'.format(data_loc,savekey) not in os.listdir(data_loc+'/plots'):
                os.mkdir('{0}/plots/{1}'.format(data_loc,savekey))

