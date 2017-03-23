import numpy as np
import pspec

class Sense(object):
    def __init__(self):
        self.Tsys = 0 #mK!
        #self.Omega_eff = 1.7 #Ali et al (in Zaki's evernote)
        self.Omega_eff = 1.74 #recalculation with capo/dcj/scripts/beam_integral.py  Omega_p^2/Omega_pp, per Appendix B of Parsons psa32
        #self.Omega_eff = 0.51**2/0.24  #effective beams if using FRF, taken from Table 1 of Parsons 2015 "beam sculpting paper"
        self.t_int = 0
        self.Npols = 2
        #note: Ndays is the effective number of days one gets after summing lst bins
        #   if the day counts vary with observed lst, then the correct number here is
        #   Ndays = 1 / sqrt(mean(1/cnt**2)) where cnt is number per lst bin vs lst
        self.Ndays = 0

        self.Nlsthours = 0
        self.Nbls = 0
        self.P_N = None
        self.Nseps = 0
        self.Nlstbins = None
        self.Nblgroups = 0
        self.Nkparfold = 2
        self.z = 0
        return
    def calc(self):
        self.X2Y = pspec.X2Y(self.z)/1e9  #515  @ z=8.4
        if self.Nlstbins is None: self.Nlstbins = self.Nlsthours*3600/self.t_int
        Nb = (self.Nbls/self.Nblgroups)
        if self.Nblgroups !=0
            self.bl_eff = Nb * np.sqrt((self.Nblgroups**2 - self.Nblgroups)/2)
        else:
            self.bl_eff = self.Nbls        
        self.bl_eff = Nb * np.sqrt((self.Nblgroups**2 - self.Nblgroups)/2)
        self.P_N = self.X2Y * self.Omega_eff * self.Tsys**2
        self.P_N /=(self.t_int * self.Ndays * self.bl_eff * self.Npols * np.sqrt(self.Nlstbins))
        self.P_N /= np.sqrt(2) #factor of 1/sqrt(2) for variance of complex variance and taking real
        self.P_N /= np.sqrt(self.Nkparfold) #fold in k//
        self.P_N /= np.sqrt(self.Nseps)
    def Delta2_N(self,k):
        if self.P_N is None: print("noise undefined until first Sense.calc()"); return 0
        return self.P_N * k**3/(2*np.pi**2)
if __name__ == "__main__":
    #an example set up for psa64 - ali et al
    S_FRF = Sense()
    S_FRF.z = 8.36
    f = pspec.z2f(S_FRF.z)*1e3

    #   Tsys
    #S_FRF.Tsys = 551e3  #set to match 21cmsense exactly
    #S.Tsys = 505e3 #Ali et al, at 164MHz
    S_FRF.Tsys = (200 + 180.*(f/180)**-2.55)*1e3 #set to match noise realization
    print "Tsys = ",S_FRF.Tsys

    S_FRF.t_int = 3414 # python ~/scripts/frf_diagnose.py -C psa6240_v003 --seps=0,1 -pxx; returns T_eff = 3414s (NEBW)
    S_FRF.Nlsthours = 8.24
    S_FRF.Ndays = 32
    #S_FRF.Nlstbins = 1  #either Nlsthours or Nlstbins must be set
    S_FRF.Nbls = 51
    S_FRF.Nseps = 3
    S_FRF.Nblgroups = 5
    S_FRF.Omega_eff = 0.51**2/0.24 #use the FRF weighted beams listed in T1 of Parsons etal beam sculpting paper
    S_FRF.calc()
    print "analytic \Delta^2(k=0.3) = ",S_FRF.Delta2_N(0.3)*2 #2sigma
