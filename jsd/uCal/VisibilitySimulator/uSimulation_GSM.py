import aipy as a
import numpy as np
import optparse, sys, os
import capo.omni as omni, capo.zsa as zsa, capo.arp as arp
import healpy as hp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy import units as u
from astropy.coordinates import SkyCoord
import scipy.constants as const

o = optparse.OptionParser()
o.set_usage('Simulate_PAPER.py [omnical .npz file to copy as a simulation] -C [calfile] -p [polarization] [other options]')
a.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('--sfreq', type='float', default=0.1,
    help='Starting frequency in GHz. Default is 0.1.')
o.add_option('--sdf', type='float', default=0.1/203,
    help='Frequency channel width in GHz. Default is 0.1/203.')
o.add_option('--nchan', type='int', default=203,
    help='Number of frequency channels. Default is 203.')
o.add_option('--NSIDE', type='int', default=256,
    help='NSIDE of GSM to use for the simulation. Default is 256.')
o.add_option('--GSMlocation', type='string', default='./GSMComponents',
    help='Location of GSM components. Default is "./GSMComponents"')
o.add_option('--nint', type='int', default=19,
    help='Number of intergrations per file. Default is 19.')
o.add_option('--inttime', type='float', default=600.0/19,
    help='Length of a single integration in seconds. Default is 600.0/19.')
o.add_option('--overwrite', action='store_true',
    help='Overwrites previous result with the same filename.')


opts,args = o.parse_args(sys.argv[1:])
if opts.cal is None: 
	opts.cal='psa6622_v003'; print 'No calfile given, defaulting to psa6622_v003'
if opts.pol is -1: 
	opts.pol = 'xx'; print 'No polarization given, defaulting to xx'
if len(args)>0: sourceFiles = args
else: sourceFiles = ['../Data/zen.2456679.35658.xx.npz']; print 'No source files given, defaulting to ../Data/zen.2456679.35658.xx.npz'
NSIDE = opts.NSIDE
NSIDE = 64

from scipy.interpolate import interp1d
class GlobalSkyModel:
    def __init__(self,GSMlocation,GSMNSIDE):
        self.NSIDE = GSMNSIDE
        self.GSMComponentsDegraded = np.asarray([np.load(GSMlocation + "/component_maps_408locked_NSIDE-" + str(GSMNSIDE) + "_Comp-" + str(comp) + ".npy") for comp in range(3)])
        self.components = np.loadtxt(GSMlocation + "/components.dat")
    def healpixMap(self, freq):
        """Returns healpix map at the given frequency (provided in MHz)."""
        temperature = np.exp(interp1d(np.log(self.components[:,0]), np.log(self.components[:,1]), kind='cubic')(np.log(freq))) #cubic spline interpolation in log(f), log(T)
        weights = np.asarray([interp1d(np.log(self.components[:,0]), self.components[:,i+2], kind='cubic')(np.log(freq)) for i in range(3)]) #cubic spline interpolation for log(f), weights
        return temperature*np.dot(weights,self.GSMComponentsDegraded)
GSM = GlobalSkyModel(opts.GSMlocation,NSIDE)


#Hardcoded!!!
blIndices = [(0, 103), (1, 4), (0, 101), (0, 62), (0, 100), (1, 13), (1, 70), (1, 56), (1, 71), (1, 59), (0, 97), (12, 43), (9, 71), (9, 59), (57, 64)]
conjugate = [(0,101), (0,62), (0,100), (0,97), (12,43), (57,64)]
aa = a.cal.get_aa(opts.cal, opts.sdf, opts.sfreq, opts.nchan)
freqs = aa.get_afreqs()
#Fixed coordinate calculations
healpixCoords = np.transpose(np.asarray(hp.pix2ang(NSIDE,np.arange(12*NSIDE**2))))
pixelRAs = healpixCoords[:,1]
pixelDecs = -healpixCoords[:,0] + np.pi/2
galCoords = SkyCoord(frame="icrs", ra=pixelRAs*u.rad, dec=pixelDecs*u.rad).transform_to("galactic")


#I should probably make this command line parameterizable
allus = np.arange(10,160,.125)


for source in sourceFiles:
    meta,_,_,_ = omni.from_npz(source)
#    times = meta['jds']
    times = [meta['jds'][0]]
    meta['SimulationOptions'] = opts.__dict__
    #vismdl = {opts.pol: {blIndex: np.zeros((len(times), opts.nchan), dtype=complex) for blIndex in blIndices}}
    visibilities = np.zeros((len(allus), len(freqs)),dtype=complex)
    
    allPixelTops, allPixelAzAlts = [],[] 
    for time in times:
        aa.set_jultime(time)
        pixelEq = a.coord.radec2eq([pixelRAs,pixelDecs])
        pixelTop = a.coord.eq2top_m(-aa.sidereal_time(),aa.lat).dot(pixelEq)
        allPixelTops.append(pixelTop)
        allPixelAzAlts.append(a.coord.top2azalt(pixelTop))

    for chan in range(opts.nchan):
        print '    Now working on channel ' + str(chan) + '...'
        GSMmapGalactic = GSM.healpixMap(freqs[chan]*1e3)
        GSMmap = hp.get_interp_val(GSMmapGalactic,-galCoords.b.radian+np.pi/2, np.asarray(galCoords.l.radian))
        convertJyToKFactor =  (a.const.c/1e2)**2 / (2 * a.const.k/1e7 * (freqs[chan]*1e9)**2 * 1e26)

        #aa.select_chans(chan)        
        
        for tIndex, time, pixelTop, pixelAzAlt in zip(range(len(times)), times, allPixelTops, allPixelAzAlts):
            #print '        JD = ' + str(time)
            aa.set_jultime(time)
            FWHM = 60.0 * np.pi/180.0 #* (.145 / freqs[chan])
            sigmaAngle = FWHM / (8*np.log(2))**.5
            sigmaTop = np.sin(sigmaAngle)
            pixelBeam = np.exp(-(pixelTop[0]**2+pixelTop[1]**2)/2/sigmaTop**2)
            pixelBeam[pixelAzAlt[1] < 0] = 0            
#            pixelBeam = pixelAzAlt[1]> 0
        
            #for blIndex in blIndices:
            for uindex, u in enumerate(allus):
                beamAndFringe = pixelBeam * np.exp(-2.0j * np.pi * np.dot([u,0.0,0.0],pixelTop))
                visibilities[uindex, chan] = np.sum(GSMmap * beamAndFringe) * 4*np.pi / len(GSMmap) / convertJyToKFactor
                
#%% Plotting!
            
plt.figure(1); plt.clf()
leg = []
for uindex, u in enumerate(allus):
    if u / 10.0 == int(u)/10: 
        plt.plot(freqs, np.real(visibilities[uindex,:] / freqs**-.55))
    #plt.plot(freqs, np.real(visibilities[uindex,:]))
        leg.append('u = '+str(u))
plt.legend(leg)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Real(Visibility) / f^(Mean SpIndex)')

#%%

for order in range(8):
#    order = 1
    errors = []
    fits = []
    for uindex, u in enumerate(allus):
        visHere =  (visibilities[uindex,:] / freqs**-.55)#freqs**np.sum(spectralIndices.dot(fluxes)/np.sum(1.0*fluxes)))
#        freqFunc = .15/freqs
        freqFunc = np.log(freqs)
        fit = np.polyfit(freqFunc, visHere, order)
        fits.append(fit)
        fitResult = np.polyval(fit, freqFunc)
        error = np.sum(np.abs(visHere - fitResult))/np.sum(np.abs(visHere))
        errors.append(error)
    print  'order:', order, " --- %.1e" % np.min(errors), "to %.1e" % np.max(errors)
#%%
fits = np.asarray(fits)
plt.figure(2); plt.clf()
for o in range(order+1):
    plt.plot(allus, np.real((-1)**o * fits[:,o])/np.mean(np.abs(fits[:,o])),'.-')
plt.legend([str(n) + ' order' for n in range(order,-1,-1)])    
#
for o in range(order+1): print np.mean(np.abs(fits[:,o]))
#
