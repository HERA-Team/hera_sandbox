#!/usr/bin/env python

#Simulation for uCal. Adapted from Carina Cheng's vis_simulaiton_v4.py and my VisibilitySimulator.py

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
o.add_option('--gaussbeam', type='float', default=0.0,
    help='Replaces the beam in the calfile with a gaussian beam in sin(alt) with given FWHM in degrees at 145 MHz. Default is 0 (use calfile beam).')
o.add_option('--prefix', type='string', default='sim',
    help='Replaces "zen" in the filename. Default "sim".')
o.add_option('--overwrite', action='store_true',
    help='Overwrites previous result with the same filename.')


opts,args = o.parse_args(sys.argv[1:])
if opts.cal is None: 
	opts.cal='psa6622_v004'; print 'No calfile given, defaulting to psa6622_v004'
if opts.pol is -1: 
	opts.pol = 'xx'; print 'No polarization given, defaulting to xx'
if len(args)>0: sourceFiles = args
else: 
    opts.overwrite = True
    sourceFiles = ['../Data/zen.2456679.35658.xx.npz']; 
    print 'No source files given, defaulting to ../Data/zen.2456679.35658.xx.npz'
NSIDE = opts.NSIDE


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

for source in sourceFiles:
    outfilename = './SimulatedVisibilities/' + source.split('/')[-1].replace('zen',opts.prefix)
    if os.path.exists(outfilename) and not opts.overwrite:
        print '\n    ' + outfilename + ' exists...skipping.\n'; continue
    omni.to_npz(outfilename, {}, {}, {}, {}) #make the file so that other programs running simultaneously skip this one
    meta,_,_,_ = omni.from_npz(source)
    times = meta['jds']
    meta['SimulationOptions'] = opts.__dict__
    vismdl = {opts.pol: {blIndex: np.zeros((len(times), opts.nchan), dtype=complex) for blIndex in blIndices}}

    #Precalcualte the topocentric and AzAlt coordinates of all pixels for all times
    allPixelTops, allPixelAzAlts = [],[] 
    for time in times:
        aa.set_jultime(time)
        pixelEq = a.coord.radec2eq([pixelRAs,pixelDecs])
        pixelTop = a.coord.eq2top_m(-aa.sidereal_time(),aa.lat).dot(pixelEq)
        allPixelTops.append(pixelTop)
        allPixelAzAlts.append(a.coord.top2azalt(pixelTop))

    #Loop over channels, regenerating the beam and GSM at each channel
    for chan in range(opts.nchan):
        print '    Now working on channel ' + str(chan) + '...'
        aa.select_chans(chan)
        GSMmapGalactic = GSM.healpixMap(freqs[chan]*1e3)
        GSMmap = hp.get_interp_val(GSMmapGalactic,-galCoords.b.radian+np.pi/2, np.asarray(galCoords.l.radian))
        convertJyToKFactor =  (a.const.c/1e2)**2 / (2 * a.const.k/1e7 * (freqs[chan]*1e9)**2 * 1e26)
        
        beamMap = a.map.Map(nside=NSIDE)
        t3 = np.array(beamMap.px2crd(px = np.arange(beamMap.npix()),ncrd=3))
        beamMap.map.map = aa[0].bm_response((t3[0],t3[1],t3[2]), pol=opts.pol[0])[0]**2
        beamMap.map.map[t3[2]<0] = 0

        #Loop over times, reinterpolating the beam at each time
        for tIndex, time, pixelTop, pixelAzAlt in zip(range(len(times)), times, allPixelTops, allPixelAzAlts):
            print '        JD = ' + str(time)
            aa.set_jultime(time)

            if opts.gaussbeam > 0:
                FWHM = opts.gaussbeam * np.pi/180.0 * (.145 / freqs[chan])
                sigmaAngle = FWHM / (8*np.log(2))**.5
                sigmaTop = np.sin(sigmaAngle)
                pixelBeam = (np.exp(-(pixelTop[0]**2+pixelTop[1]**2)/2/sigmaTop**2)) * pixelTop[2] > 0
                #pixelBeam = pixelTop[2] > 0
            else:
                pixelBeam = hp.get_interp_val(beamMap.map.map, np.pi/2-pixelAzAlt[1], pixelAzAlt[0])

            #Loop over baselines
            for blIndex in blIndices:
                bl = np.round(aa.get_baseline(blIndex[0],blIndex[1]) * a.const.c / 1.5e12) / a.const.c * 1.5e12
                beamAndFringe = pixelBeam * np.exp(-2.0j*np.pi*freqs[chan] * np.dot(bl,pixelTop))
                vismdl[opts.pol][blIndex][tIndex,chan] = np.sum(GSMmap*beamAndFringe) * 4*np.pi / beamMap.npix() / convertJyToKFactor

    omni.to_npz(outfilename, meta, {}, vismdl, {})
    print '\nResults written to ' + outfilename + '\n'







