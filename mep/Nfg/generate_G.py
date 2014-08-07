#! /usr/bin/env python

import aipy as ap, numpy as np, sys
import capo as C

def get_bls(calfile,na = 32,freqs = np.arange(.1,.2,.01)):
    aa = ap.cal.get_aa(calfile, np.array([.150])) #get antenna array
    baselines = np.zeros([(na*na-na)/2,3])
    ll=0
    for ii in np.arange(na):
        for jj in np.arange(ii):
            bx,by,bz = aa.get_baseline(ii,jj,'z') 
            #the baseline for antennas ii and jj 
            #print ii,jj,[bx,by,bz]
            baselines[ll] = [bx,by,bz]
            ll+=1
    return baselines

def gauss(sig,x,x0=0.,y=0.,y0=0.):
    return 1/(2*np.pi*sig*sig)*np.exp(-((x-x0)**2)/(2*sig*sig))*np.exp(-((y-y0)**2)/(2*sig*sig))

if __name__=='__main__': 
    calfile = sys.argv[1]
    na = int(sys.argv[2])
    nside = int(sys.argv[3])
    npix = 12 * nside * nside
    freq = float(sys.argv[4]) / 1000.
    beamSig = float(sys.argv[5])
    outputFname = sys.argv[6]

    baselines = get_bls(calfile,na=na,freqs=np.array([freq,]))
    numBl = len(baselines)
    dOmega = 4 * np.pi / npix

    beamMap = ap.map.Map(nside)
    directionVects_xyz = beamMap.map.px2crd(np.array([i for i in range(npix)]))
    directionVects_xyz = np.array(directionVects_xyz).T
    directionVects_thetas = beamMap.map.px2crd(np.array([i for i in range(npix)]),ncrd=2)[0]
    
    primaryBeam = gauss(beamSig,directionVects_thetas)

    Gmatrix = np.zeros((numBl,npix),dtype=complex)
    for j,(nhat,beamVal) in enumerate(zip(directionVects_xyz,primaryBeam)):
        for i,bl in enumerate(baselines):
            Gmatrix[i,j] = np.exp(-2j*np.pi*freq*np.dot(nhat,bl))
        Gmatrix[:,j] *= beamVal
    Gmatrix *= dOmega
    
    np.save(outputFname,Gmatrix)