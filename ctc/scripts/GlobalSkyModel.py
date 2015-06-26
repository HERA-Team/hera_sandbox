# SUPPORTING CLASS FOR JOINT MAPMAKING AND POWER SPECTRUM PIPELINE
# original script by Josh Dillon
# edited to be Aipy-compatible by Carina Cheng

import numpy as np
#import healpy as hp
import aipy
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class GlobalSkyModel:

        #Takes the frequency (in MHz), the location of the HEALPIX .fits files (which are in ~/capo/ctc/images/gsm/components/) and the HealPIX NSIDE desired (powers of 2, 8 through 512, are available)
    def __init__(self,freq,GSMlocation,GSMNSIDE):
        #Compute GSM from 3 principal components appropriately weighted        
        self.freq = freq
        self.NSIDE = GSMNSIDE        
        #GSMComponents = np.asarray([hp.read_map(GSMlocation + "gsm" + str(i+1) + ".fits" + str(GSMNSIDE),verbose=False) for i in range(3)])
        GSMComponents = np.asarray([aipy.map.Map(fromfits = GSMlocation + "gsm" + str(i+1) + ".fits" + str(GSMNSIDE)).map.map for i in range(3)])
        components = np.loadtxt(GSMlocation + "components.dat")
        temperature = 10**(interp1d(np.log10(components[:,0]), np.log10(components[:,1]), kind='cubic')(np.log10(freq))) #cubic spline in log(T)
        weights = np.asarray([interp1d(components[:,0], components[:,i+2], kind='linear')(freq) for i in range(3)]) #linear interpolation for weights
        self.map = temperature*np.dot(weights,GSMComponents)
