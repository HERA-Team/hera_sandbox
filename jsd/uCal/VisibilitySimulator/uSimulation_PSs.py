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


allus = np.arange(10,160,.125)


#allus = np.arange(100,102,.1)

#pixelDecs = np.asarray([-25])/180.0*np.pi
#pixelRAs = np.asarray([90])/180.0*np.pi
#fluxes = np.asarray([100])
pixelDecs = np.asarray([-30, -35, -20, -10, 5, -30, -50, -40, -60, -20])/180.0*np.pi
pixelRAs = np.asarray([90, 60, 120, 76, 91, 44, 50, 70, 30, 10])/180.0*np.pi
fluxes = np.asarray([100, 49, 125, 60, 20, 99, 22, 102, 5, 30])*1.0
#spectralIndices = np.ones(len(fluxes)) * -2.55
spectralIndices = np.asarray([-2.45, -2.3, -2.7, -2.4, -2.55, -2.63, -2.51, -2.3, -2.9, -2.5]) + 2

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
        #aa.select_chans(chan)
        
        for tIndex, time, pixelTop, pixelAzAlt in zip(range(len(times)), times, allPixelTops, allPixelAzAlts):
            # print '        JD = ' + str(time)
            aa.set_jultime(time)
            # pixelBeam = hp.get_interp_val(beamMap.map.map, np.pi/2-pixelAzAlt[1], pixelAzAlt[0])
            FWHM = 60.0 * np.pi/180.0 * (.145 / freqs[chan])
            sigmaAngle = FWHM / (8*np.log(2))**.5
            sigmaTop = np.sin(sigmaAngle)
            #pixelBeam = np.exp(-(pixelTop[0]**2+pixelTop[1]**2)/2/sigmaTop**2) * (pixelAzAlt[1]> 0)
            pixelBeam = pixelAzAlt[1]> 0
        
            #for blIndex in blIndices:
            for uindex, u in enumerate(allus):
                #bl = np.round(aa.get_baseline(blIndex[0],blIndex[1]) * a.const.c / 1.5e12) / a.const.c * 1.5e12
                beamAndFringe = pixelBeam * np.exp(-2.0j * np.pi * np.dot([u,0.0,0.0],pixelTop))
                #beamAndFringe = pixelBeam * np.exp(-2.0j*np.pi* freqs[chan] * np.dot(bl,pixelTop))
                visibilities[uindex, chan] = np.sum(fluxes*beamAndFringe * freqs[chan]**spectralIndices) #* 4*np.pi / beamMap.npix() / convertJyToKFactor
                #vismdl[opts.pol][blIndex][tIndex,chan] = visibility
                
#%% Plotting!
            
plt.figure(1); plt.clf()
for uindex, u in enumerate(allus):
    if u < 12: plt.plot(freqs, np.real(visibilities[uindex,:] / freqs**np.sum(spectralIndices.dot(fluxes)/np.sum(fluxes))))
    #plt.plot(freqs, np.real(visibilities[uindex,:]))
plt.legend(['u = '+str(u) for u in allus])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Real(Visibility) / f^(Mean SpIndex)')

#%%

for order in range(8):
#    order = 1
    errors = []
    fits = []
    for uindex, u in enumerate(allus):
        visHere =  (visibilities[uindex,:] / freqs**np.sum(spectralIndices.dot(fluxes)/np.sum(1.0*fluxes)))
#        freqFunc = .15/freqs
        freqFunc = np.log(freqs/.145)
        fit = np.polyfit(freqFunc, visHere, order)
        fits.append(fit)
        fitResult = np.polyval(fit, freqFunc)
        error = np.sum(np.abs(visHere - fitResult))/np.sum(np.abs(visHere))
        errors.append(error)
    print order, " :  %.1e" % np.min(errors), "to %.1e" % np.max(errors)
    
plt.figure(1000); plt.clf()
plt.plot(np.abs(fits))
##%%
fits = np.asarray(fits)
plt.figure(2); plt.clf()
for o in range(order+1):
    plt.plot(allus, np.real((-1)**o * fits[:,o])/np.mean(np.abs(fits[:,o])),'.-')
plt.legend([str(n) + ' order' for n in range(order,-1,-1)])    
#
for o in range(order+1): print np.mean(np.abs(fits[:,o]))
#
