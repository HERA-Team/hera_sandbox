"""
vis_sim.py
=============

Routines for visibility simulation using the
Global Sky Model (GSM) and an antenna beam model
"""

# Load Modules
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import healpy as hp
import os
import sys
from scipy import interpolate
import astropy.io.fits as fits
from astropy.time import Time
import ephem
from hera_cal import omni
import aipy
from collections import OrderedDict
from astropy.stats import biweight_midvariance
import copy
from pyuvdata import UVData, UVBeam
import itertools


class Helper(object):
    """
    Helper
    -------

    Helper routines
    """
    def print_message(self, message, type=0, verbose=True):
        """
        print message to stdout
        """
        if verbose == True:
            if type == 0:
                print("\n%s" % message)
            elif type == 1:
                print("\n%s\n%s" % (message, '-'*40))

    def JD2LST(self, JD, longitude=None, loc=None):
        """
        use astropy to convert from JD to LST
        provide either:
        loc : ephem Observer instance
        longitude : float, longitude in degrees East
        """
        if loc is not None:
            t = Time(JD, format="jd")
            lst = t.sidereal_time('apparent', longitude=loc.lon*180/np.pi)
        elif longitude is not None:
            t = Time(JD, format="jd")
            lst = t.sidereal_time('apparent', longitude=longitude)
        return lst.value

    def healpix_interp(self, maps, map_nside, theta, phi, nest=False):
        """
        healpix map bi-linear interpolation

        Input:
        maps : ndarray, shape=(Nmaps, Npix)
            array containing healpix map(s)

        map_nside : int
            nside parameter for maps

        theta : ndarray, shape=(Npix)
            1D array containing co-latitude of
            desired inteprolation points in radians

        phi : ndarray, shape=(Npix)
            1D array containing longitude of
            desired interpolatin points in radians

        """
        if nest == True:
            r = hp._healpy_pixel_lib._get_interpol_nest(map_nside, theta, phi)
        else:
            r = hp._healpy_pixel_lib._get_interpol_ring(map_nside, theta, phi)
        p=np.array(r[0:4])
        w=np.array(r[4:8])
        if maps.ndim == 2:
            return np.einsum("ijk,jk->ik", maps[:, p], w)
        elif maps.ndim == 3:
            return np.einsum("hijk,jk->hik", maps[:, :, p], w)
        elif maps.ndim == 4:
            return np.einsum("ghijk,jk->ghik", maps[:, :, :, p], w)

    def rotate_map(self, nside, rot=None, coord=None, theta=None, phi=None, interp=False,
                   inv=True):
        """
        rotate healpix map between coordinates and/or in longitude and latitude

        nside : int, nside resolution of map

        rot : list, length=2
            rot[0] = longitude rotation (radians)
            rot[1] = latitude rotation (radians)

        coord : list, length=2
            transformation between coordinate systems
            see healpy.Rotator for convention

        theta : co-latitude map in alt-az coordinates

        phi : longitude map in alt-az coordinates

        interp : bool, default=False
            if True, use interpolation method (healpy.get_interp_val)
            if False, use slicing method (healpy.ang2pix)

        inv : bool, default=True
            keyword to feed hp.Rotator object
            to go from galactic to topocentric, inv=True
            to go from topocentric to galactic, inv=False
        """
        # if theta and phi arrays are not fed, build them
        if theta is None or phi is None:
            theta, phi = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

        # get rotation
        rot_theta, rot_phi = hp.Rotator(rot=rot, coord=coord, inv=inv, deg=False)(theta, phi)

        if interp is False:
            # generate pixel indices array
            pix = np.arange(hp.nside2npix(nside))[hp.ang2pix(nside, rot_theta, rot_phi)]
            return pix

        else:
            return rot_theta, rot_phi


class GlobalSkyModel(Helper):
    """
    GlobalSkyModel
    --------------

    class for handling of global sky model data
    """
    def __init__(self, gsmfile, sky_nside=None, freqs=None, verbose=False, onepol=False):
        """
        load global sky model
        and optionally interpolate to specific nside
        spatial resolution and frequency resolution

        Input:
        -------
        gsmfile : str
            path to gsm data file as a .npz file
            with "sky" as multifrequency sky and
            "freqs" as frequency channels in MHz.
            sky is assumed to be in galactic coordinates.

        sky_nside : int, default=None
            downsample GSM to desired healpix nside
            resolution (not recommended).
            If None, default is what the gsmfile provides

        freqs : ndarray, default=None
            frequency channels of sky model in MHz. If None, default
            is what the gsmfile provides

        onepol : bool, default=False
            if True, keep only stokes I map
            if False, keep full Jones matrix

        Result:
        -------
        self.sky_models : ndarray, shape=(2, 2, Nfreqs, Npix)
            first two axes are for 4-pol stokes parameters
            [ [I+Q, U+iV], [U-iV, I-Q] ]
            (multi-frequency) global sky model in healpix
            using RING ordering in units of MJy-str

        self.freqs : ndarray, shape=(Nfreqs,)
            frequency of sky maps in MHz

        self.theta : ndarray, shape=(Npix,)
            angle of latitude of healpix map in radians
            in galactic coordinates

        self.phi : ndarray, shape=(Npix,)
            angle of longitude of healpix map in radians
            in galactic coordinates
        """

        # Load GSM Sky 
        GSM_data    = np.load(gsmfile)
        sky_models  = GSM_data['sky']
        sky_freqs   = GSM_data['freqs']
        sky_freq_ax = 2

        # check sky_models are appropriate shape
        if len(sky_models.shape) != 4:
            raise ValueError("sky_models.shape = {} when it should have len = 4".format(sky_models.shape))

        # check for onepol
        if onepol is True:
            sky_models = (np.trace(sky_models, axis1=0, axis2=1) / 2.0)[np.newaxis]
            sky_freq_ax -= 1

        # 1D interpolation of frequency axis
        if freqs is not None:
            sky_models = interpolate.interp1d(sky_freqs, sky_models, axis=sky_freq_ax)(freqs)
            sky_freqs = freqs

        # get theta and phi arrays for default data
        default_sky_nside = hp.npix2nside(sky_models.shape[sky_freq_ax+1])
        theta, phi = hp.pix2ang(default_sky_nside, np.arange(12*default_sky_nside**2), lonlat=False)

        # spatially interpolate healpix if desired
        if sky_nside is not None:
            # down sample
            theta, phi = hp.pix2ang(sky_nside, np.arange(12*sky_nside**2))
            sky_models = self.healpix_interp(sky_models, default_sky_nside, theta, phi)
            default_sky_nside = sky_nside

        # Assign variables to class
        self.sky_nside  = default_sky_nside
        self.sky_npix   = 12 * default_sky_nside**2
        self.sky_models = sky_models
        self.sky_freqs  = sky_freqs
        self.sky_theta  = theta
        self.sky_phi    = phi
        self.onepol     = onepol
        self.sky_freq_ax = sky_freq_ax

    def GSM_freq_interp(self, freqs):
    	"""
    	1D interpolation of sky models along frqeuency axis

    	freqs : ndarray
    		frequency channels in MHz
    	"""
    	self.sky_models = interpolate.interp1d(self.sky_freqs, self.sky_models, axis=self.sky_freq_ax)(freqs)
    	self.sky_freqs = freqs

    def plot_sky(self, loc, skymap, ax=None, log10=True, res=300, axoff=True, cbar=True,
    				save=False, fname=None, plot_kwargs={'cmap':'viridis'},
    				cbar_kwargs={}, basemap=True, rot=None):
        """
		Plot Global Sky Model in orthographic coordinates given observer
		location and date

		Input:
		------
		loc : ephem location object
			relevant subattributes are loc.lon, loc.lat and loc.date

		skymap : ndarray, shape=(Npix,)
			GSM map in healpix RING ordered

		ax : matplotlib axis object, default=None
			feed a previously defined axis if desired
			else create a new figure and axis object

		log10 : bool, default=True
			take the log10 of the map before plotting

		res : int, default=300
			polar coordinates pixel resolution
			plotting gets slow when res > 500

		axoff : bool, default=True
			turn off axes tick labels

		cbar : bool, default=True
			make a colorbar

		save : bool, default=False
			save image to file

		fname : str, default=None
			filename of image to save

		plot_kwargs : dictionary
			keyword arguments to feed plotting routine

		cbar_kwargs : dictionary
			kwargs to feed colorbar routine
		
		Output:
		-------
		if ax is None:
			outputs matplotlib.pyplot.figure object
        """
        # rotate map
        if rot is not None:
            rot = self.rotate_map(self.sky_nside, rot=rot)
            skymap = skymap[rot]

        # get ra dec of location
        obs_ra, obs_dec = loc.radec_of(0, np.pi/2.0)

        # get rotation sorting array
        rot = self.rotate_map(self.sky_nside, rot=[obs_ra, obs_dec-np.pi/2], coord=['G', 'C'], theta=self.sky_theta, phi=self.sky_phi)

        # get rotated sky
        rot_sky = skymap[rot]

        # Get polar theta and r
        omega, r = np.meshgrid(np.linspace(0,2*np.pi,res), np.linspace(0,np.pi/2.0,res))

        # sample sky at polar coordinates
        obs_sky_polar = hp.get_interp_val(rot_sky, r.ravel(), omega.ravel()).reshape(res,res)

        # rotate omega
        omega  = np.unwrap(omega + np.pi/2.0)

        # mirror about x-axis
        i = r * np.cos(omega) + 1j * r * np.sin(omega)
        i = np.conj(i)
        omega = np.unwrap(np.angle(i))
        r = np.abs(r)

        # create custom fig, ax if necessary
        custom_fig = False
        if ax is None:
            custom_fig = True
            if basemap == True:
                fig, ax = plt.subplots(1, figsize=(6,6))
            else:
                fig, ax = plt.subplots(1, figsize=(6,6), subplot_kw=dict(projection='polar'))

        # log data if desired
        if log10 == True:
            obs_sky_polar = np.log10(obs_sky_polar)

        # use basemap
        if basemap == True:
            bmap = Basemap(projection='ortho',lat_0=90,lon_0=-90, ax=ax)

            # turn r from co-latitude to latitude
            r = np.abs(r-np.pi/2)

            # get x, y arrays for basemap
            x, y = bmap(omega*180/np.pi, r*180/np.pi)

            # plot
            cax = bmap.pcolor(x, y, obs_sky_polar, **plot_kwargs)

        # use polar
        else:
            # plot
            cax = ax.pcolormesh(omega, r, obs_sky_polar, **plot_kwargs)

        # turn axis off
        if axoff == True:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

        # add colorbar
        if cbar == True:
            plt.colorbar(cax, orientation='vertical', **cbar_kwargs)

        if save == True:
            plt.savefig(fname, dpi=200, bbox_inches='tight')

        if custom_fig == True:
            return fig


class Beam_Model(Helper):
    """
    Beam_Model
    ----------

    class for handling of beam models
    """
    def __init__(self, beamfile, loc=None, freqs=None, pool=None, verbose=False,
                 mask_hor=True, beam_nside=None, onepol=False, pol=0, beam_type='intensity_sq'):
        """
        Load and configure beam models

        Input:
        ------
        beamfile : str
            pyuvdata beamfits file, holding a healpix beam in topocentric coordinates

        freqs : ndarray, default=None
            frequency channels in MHz, defaults to
            channels in steps of 1 MHz

        beam_data_path : str
        	path to beam data in fits format with 2 HDUs
        	with healpix beam models in 0th hdu w/ shape=(Npix, Nfreq)
        	and beam frequencies in 1st hdu w/ shape=(Nfreq,)

        beam_nside : int
            interpolate beam onto beam_nside

        onepol : bool, default=False
            if True: use only one feed polarization (zeroth axis)

        pol : int, default=0, option = 0 or 1
            if onepol is True, which pol index to use in beam_models 

        beam_type : str, default='power', options=['intensity', 'intensity_sq']
            specifies the beam_type of the map. the measurement equation
            reads V = A1.T * I * A2, where A is the beam intensity and I is the sky
            brightness. If we assume A1 == A2 and A3 == A2**2, then this
            simplifies to V = A3 * I, where A1 and A2 have type of 'intensity'
            and A3 has type of 'intensity_sq', i.e. intensity_sq = intensity**2.
            this simulation needs a beam type of 'intensiy', so if fed
            'intensity_sq', it will take the square root of the maps.

        Result:
        -------
        self.beam_models : masked array, shape=(Nfreq, Npix)
            beam model map in healpix with sub-horizon masked out

        self.beam_freqs : ndarray, shape=(Nfreqs,)
            ndarray containing frequencies of beam in MHz

        self.beam_theta : ndarray, shape=(Npix,)
            ndarray containing co-latitude in radians of
            healpix map

        self.beam_phi : ndarray, shape=(Npix,)
            ndarray containing longitude in radians of
            healpix map
        """
        # Assign telescope coordinates and location
        if loc is None:
            loc = ephem.Observer()
        self.loc = loc

        # Load Models in healpix from fits file
        uvb         = UVBeam()
        uvb.read_beamfits(beamfile)
        beam_models = uvb.data_array[0,0]
        beam_freqs  = uvb.freq_array.squeeze() / 1e6
        pol_arr     = uvb.polarization_array

        # ensure beam type is intensity
        if beam_type == 'intensity_sq':
            beam_models = beam_models ** (1./2)

        # ensure beam_models.shape
        if len(beam_models.shape) < 2 or len(beam_models.shape) > 3:
            raise ValueError("beam_models.shape is {}, but must have len of either 2 or 3".format(beam_models.shape))

        # select onepol
        if onepol is True:
            beam_models = beam_models[pol][np.newaxis]
            pol_arr = np.array(pol_arr[pol])[np.newaxis]
        else:
            if beam_models.shape[0] != 2:
                raise ValueError("beam_models.shape[0] must be 2, but is {}".format(beam_models.shape[0]))

        # 1D interpolation of frequency axis if desired
        if freqs is not None:
            beam_models = interpolate.interp1d(beam_freqs, beam_models, axis=1,
                                               bounds_error=False, fill_value='extrapolate')(freqs)
            beam_freqs = freqs

        # Get theta and phi arrays
        default_beam_nside = uvb.nside
        default_beam_npix = uvb.Npixels
        beam_theta, beam_phi = hp.pix2ang(default_beam_nside, np.arange(default_beam_npix), lonlat=False)

        # spatially interpolate healpix if desired
        if beam_nside is not None:
            # down sample
            default_beam_nside = beam_nside
            default_beam_npix = hp.nside2npix(default_beam_nside)
            beam_theta, beam_phi = hp.pix2ang(default_beam_nside, np.arange(default_beam_npix))
            beam_models = np.array(map(lambda x: map(lambda y: hp.get_interp_val(y, beam_theta, beam_phi), x), beam_models))

        # make sure boresight is normalized to 1
        beam_models /= np.max(beam_models, axis=-1)[:, :, np.newaxis]

        # mask beam models below horizon
        if mask_hor is True:
            mask = (beam_theta > np.pi/2)
            beam_models[:, :, mask] = 0.0

        # assign vars to class
        self.beam_nside  = default_beam_nside
        self.beam_npix   = default_beam_npix
        self.beam_models = beam_models
        self.beam_freqs  = beam_freqs
        self.beam_theta  = beam_theta
        self.beam_phi    = beam_phi
        self.beam_pols   = pol_arr

    def beam_freq_interp(self, freqs):
        """
        1D interpolation of beam models along frqeuency axis

        freqs : ndarray
            frequency channels in MHz
        """
        self.beam_models = interpolate.interp1d(self.beam_freqs, self.beam_models, axis=1)(freqs)
        self.beam_freqs = freqs

    def project_beams(self, JD, sky_theta, sky_phi, beam_models=None, obs_lat=None, obs_lon=None,
    					freqs=None, output=False, pool=None, interp=False):
        """
        Project beam models into healpix galactic coordinates
        given observer location, observation date and sky models
        and interpolate onto sky model healpix nside resolution

        Input:
        ------
        JD : float
            Julian date of observation in J2000 epoch

        sky_theta : ndarray, shape=(Npix,)
            co-latitude in radians of sky healpix map
            in galactic coordinates

        sky_phi : ndarray, shape=(Npix,)
            longitude in radians of sky healpix map
            in galactic coordinates

		beam_models : list, shape=(Npol, Nfreq, Npix)
			set of beam models, default=None
			in which case it will use self.beam_models

        obs_lat : str or float
            observer's latitude on Earth
            if str type: format "deg:hour:min"
            if float type: format is radians

        obs_lon : str or float
            observer's longitude on Earth
            if str type: format "deg:hour:min"
            if float type: format is radians

        Result:
        -------
        self.sky_beam_models : masked array, shape=(Nfreq, Npix)
            beam models projected onto galactic coordaintes
            and interpolated onto sky model healpix resolution
            at each frequency with sub-horizon masked out

        self.s_hat : ndarray, shape=(Npix,3)
        	ndarray containing pointing unit vector in observer's
        	cartesian frame for each healpix pixel of beam
        """
        # get beam models
        if beam_models is None:
        	beam_models = self.beam_models 

        # Assign coordinates and date of observation
        self.loc.date = Time(JD, format='jd').datetime
        if obs_lon is not None:
            self.loc.lon = obs_lon
        if obs_lat is not None:
            self.loc.lat = obs_lat

        # get ra/dec of zenith
        obs_ra, obs_dec = self.loc.radec_of(0, np.pi/2.0)

        # get rotation sorting array
        if interp is True:
            rot_theta, rot_phi = self.rotate_map(self.sky_nside, rot=[obs_ra, obs_dec-np.pi/2], coord=['G', 'C'], theta=sky_theta, phi=sky_phi, inv=False, interp=interp)
            self.sky_beam_models = self.healpix_interp(beam_models, self.beam_nside, rot_theta, rot_phi)

        else:
            rot = self.rotate_map(self.sky_nside, rot=[obs_ra, obs_dec-np.pi/2], coord=['G', 'C'], theta=sky_theta, phi=sky_phi, inv=False, interp=interp)
            self.sky_beam_models = beam_models[:, :, rot]

        self.sky_beam_freqs = self.beam_freqs

        # Interpolate frequency axis if desired
        if freqs is not None:
        	self.sky_beam_models = map(lambda sbm: interpolate.interp1d(self.beam_freqs, sbm, axis=0)(freqs), self.sky_beam_models)
        	self.sky_beam_freqs = freqs

        if output == True:
        	return self.sky_beam_models


    def plot_beam(self, beam, ax=None, log10=False, dBi=False, res=300, axoff=True, cbar=False,
    				contour=True, levels=None, nlevels=None, label_cont=False,
    				save=False, fname=None, basemap=True, verbose=False,
    				plot_kwargs={'cmap':'YlGnBu_r','linewidths':0.75,'alpha':0.8},
    				cbar_kwargs={}, rot=[0, np.pi/2]):
        """
		Plot beam model in orthographic coordinates

		Input:
		------
		beam : ndarray, shape=(Npix,)
			beam response in healpix RING ordered

		ax : matplotlib axis object, default=None
			feed a previously defined axis if desired
			else create a new figure and axis object

		log10 : bool, default=False
			take the log10 of the map before plotting

        dBi : bool, default=False
            express data in deciBels of intensity
            in other words: log10(data) * 10

		res : int, default=300
			polar coordinates pixel resolution
			Beware: plotting gets slow when res > 500

		axoff : bool, default=True
			turn off axes tick labels

		cbar : bool, default=True
			make a colorbar

		contour : bool, default=True
			if True: contour plot
			else: pcolormesh plot

		levels : list, default=None
			levels for contours

		nlevels : int, default=None
			if levels is None, try to 
			make custom levels with nlevels
			if None and levels is None, contour()
			will make its own levels by default

		label_cont : bool, default=False
			label contours with level

		save : bool, default=False
			save image to file

		fname : str, default=None
			filename of image to save

		plot_kwargs : dictionary
			keyword arguments to feed plotting routine

		cbar_kwargs : dictionary
			kwargs to feed colorbar routine
		
		Output:
		-------
		if ax is None:
			outputs matplotlib.pyplot.figure object
        """
        if dBi == True and log10 == True:
            log10 = False
            if verbose == True:
                print "...both log10 dBi is True, using dBi"

        # rotate beam
        if rot is not None:
            rot = self.rotate_map(self.beam_nside, rot=rot)
            beam = beam[rot]

        # get beam_theta (co-latitude) and beam_phi (longitude)
        nside = int(np.sqrt(len(beam)/12.0))
        beam_theta, beam_phi = hp.pix2ang(nside, np.arange(len(beam)), lonlat=False)

        # rotate observed beam upwards
        theta_pol, phi_pol = hp.Rotator(rot=[0,np.pi/2,0],deg=False,inv=False)(beam_theta,beam_phi)
        rot_beam = hp.get_interp_val(beam, theta_pol, phi_pol)

        # Get polar theta and r
        omega, r = np.meshgrid(np.linspace(0,2*np.pi,res), np.linspace(0,np.pi/2,res))

        # sample sky at polar coordinates
        beam_polar = hp.get_interp_val(rot_beam, r.ravel(), omega.ravel()).reshape(res,res)

        # rotate omega by np.pi/2
        omega  = np.unwrap(omega + np.pi/2.0)

        # mirror about x-axis
        i = r * np.cos(omega) + 1j * r * np.sin(omega)
        i = np.conj(i)
        omega = np.unwrap(np.angle(i))
        r = np.abs(r)

        # create custom fig, ax if necessary
        custom_fig = False
        if ax is None:
            custom_fig = True
            if basemap == True:
                fig, ax = plt.subplots(1, figsize=(6,6))
            else:
                fig, ax = plt.subplots(1, figsize=(6,6), subplot_kw=dict(projection='polar'))

        # log data if desired
        if log10 == True:
            beam_polar = np.log10(beam_polar)

        if dBi == True:
            beam_polar = np.log10(beam_polar) * 10

        # use basemap
        if basemap == True:
            bmap = Basemap(projection='ortho',lat_0=90,lon_0=-90, ax=ax)

            # turn r from co-latitude to latitude
            r = np.abs(r-np.pi/2)

            # get x, y arrays for basemap
            x, y = bmap(omega*180/np.pi, r*180/np.pi)

            # plot
            if contour == True:
                # get contour levels
                if levels is None and nlevels is not None:
                    maxval = beam_polar.max() * 0.95
                    levels = [maxval/i for i in np.arange(1.0, nlevels+1)][::-1]
                cax = bmap.contour(x, y, beam_polar, levels=levels, **plot_kwargs)
            else:
                cax = bmap.pcolor(x, y, beam_polar, **plot_kwargs)
        # use polar
        else:
            # plot
            if contour == True:
                # get contour levels
                if levels is None and nlevels is not None:
                    maxval = beam_polar.max() * 0.95
                    levels = [maxval/i for i in np.arange(1.0, nlevels+1)][::-1]
                cax = ax.contour(omega, r, beam_polar, levels=levels, **plot_kwargs)
            else:
                cax = ax.pcolormesh(omega, r, beam_polar, **plot_kwargs)

        # turn axis off
        if axoff == True:
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)

        # add colorbar
        if cbar == True:
            plt.colorbar(cax, orientation='horizontal', **cbar_kwargs)

        # label contours
        if label_cont == True:
            ax.clabel(cax, inline=1, fontsize=10)

        if save == True:
            plt.savefig(fname, dpi=200, bbox_inches='tight')

        if custom_fig == True:
            return fig

class Vis_Sim(GlobalSkyModel, Beam_Model):
    """
    Vis_Sim
    -------

    class for handling visibility simulations

    """

    def init_GSM(self, gsmfile, **kwargs):
    	"""
    	initialize GSM. See GlobalSkyModel class doc-string for info on arguments
    	"""
    	super(Vis_Sim, self).__init__(gsmfile, **kwargs)

    def init_beam(self, beamfile, **kwargs):
        """
        initialize beam. See Beam_Model class doc-string for info on arguments
        """
        # initialize Beam_Model class
        super(GlobalSkyModel, self).__init__(beamfile, **kwargs)

        # get polarization parameters
        self.pols = self.beam_pols.copy()
        self.Npols = len(self.pols)

        # set fiducial beam models
        self.fid_beam_models = copy.deepcopy(self.beam_models)

        # Initialize principal components
        self.Npcomps = 1
        self.pcomps  = np.ones((self.Npols, self.Nfreqs, self.Npcomps, self.beam_npix))

        # Initialize beam coefficients for each antenna
        self.beam_coeffs = np.zeros((self.red_info.nAntenna, self.Npols, self.Nfreqs, self.Npcomps))


    def __init__(self, calfile, gsmfile, beamfile, info_kwargs={}, sky_nside=None, freqs=None,
                 onepol=False, pol=0, verbose=False):
        """
		Redundant Array Visibility Simulation with the Global Sky Model
		Inherits from classes <GlobalSkyModel> and <Beam_Model>

		Input:
		------
		calfile : str
			name of an aipy calfile w/o .py suffix
			within your current path

        gsmfile : str
            data for gsmfile in .npz format

        beamfile : str
            data for beamfile as a pyuvdata .beamfits file

        info_kwargs : dictionary, default={}
            kwargs to feed omni.aa_to_info() function

        sky_nside : int, default=None
            nside healpix resolution parameter for sky model

        freqs : ndarray, default=None
            desired visibility frequency channels in MHz
            by default uses sky model frequency channels

        onepol : bool, default=False
            if True: solve only one auto-pol (faster)
                takes trace/2 of sky models (stokes I) and [pol] of beam models
            if False: solve full jones matrix for 4-pol terms (slower)

        pol : int, default=0
            if onepol is True, which polarization index to use for beam models (0 or 1)

        verbose : bool, default=False
            print out progress

        """
        self.print_message("...initializing visibility simulation", type=1, verbose=verbose)

        # Initialize GSM
        self.init_GSM(gsmfile, sky_nside=sky_nside, freqs=freqs, verbose=verbose, onepol=onepol)

        # Assign variables
        self.calfile     = calfile
        self.info_kwargs = info_kwargs
        if freqs is None:
            self.freqs   = self.sky_freqs
        else:
            self.freqs   = freqs
        self.Nfreqs      = len(self.freqs)

        # Get Redundant Info
        self.print_message("...getting redundant info", type=1, verbose=verbose)
        self.AntArr   = aipy.cal.get_aa(calfile, self.freqs)
        self.red_info = omni.aa_to_info(self.AntArr, **info_kwargs)
        self.red_bls  = np.array(self.red_info.get_reds()) 
        self.ant_nums = np.unique(np.concatenate(self.red_bls.tolist()).ravel())
        self.Nants    = len(self.ant_nums)
        self.ant_pos  = OrderedDict([(str(i), self.red_info.antloc[self.red_info.ant_index(i)]) for i in self.ant_nums])
        self.ant2ind  = OrderedDict(zip(np.array(self.ant_nums,str), np.arange(self.Nants)))
 
        # Assign location info
        self.loc           = ephem.Observer()
        self.loc.lon       = self.AntArr.lon
        self.loc.lat       = self.AntArr.lat
        self.loc.elev      = self.AntArr.elev

        # Initialize fiducial beam models if beamfile
        self.print_message("...initializing beam models", type=1, verbose=verbose)
        self.init_beam(beamfile, loc=self.loc, pol=pol, onepol=onepol, freqs=self.freqs,
                       verbose=verbose, beam_nside=sky_nside)

    def build_beams(self, ant_inds=None, output=False, one_beam_model=True):
    	""""
    	build each antenna's beam model from fiducial model and principal components

    	Input:
    	------
    	output : bool, default=False
    		if True, return results

    	Result:
    	-------
    	self.ant_beam_models : dictionary
    		contains each antenna name in str as key
    		holding healpix beam model
    	"""
        # get antenna indices
        if ant_inds is None:
            ant_inds = np.arange(self.Nants)            

        # get each antennas beam model
        if one_beam_model == True:
            self.ant_beam_models = self.fid_beam_models.copy()
        else:
            self.ant_beam_models = np.array(map(lambda i: np.einsum('ijkl,ijk->ijl', self.pcomps, self.beam_coeffs[i]) + self.fid_beam_models, ant_inds))

    	if output == True:
    		return self.ant_beam_models

    def sim_obs(self, bl_array, JD_array, pool=None,
                write_miriad=False, fname=None, clobber=False, one_beam_model=True,
                interp=True):
        """
		Simulate a visibility observation of the sky

		Input:
		------
		bl_array : list, shape=(Nbls, 2), entry_format=tuple(int, int)
		    list of baselines (antenna pairs)
		    in (ant1, ant2) with type(ant1) == int

		JD_array : list, shape=(Ntimes,), dtype=float
			list of Julian Dates for observations
    	"""
        # get array info
        Nbls         = len(bl_array)
        str_bls      = [(str(x[0]), str(x[1])) for x in bl_array]
        Ntimes       = len(JD_array)
        rel_ants     = np.unique(np.concatenate(bl_array))         # relevant antennas
        rel_ant2ind  = OrderedDict(zip(np.array(rel_ants,str), np.arange(len(rel_ants))))

        # get antenna indices
        ant_inds = map(lambda x: self.ant2ind[str(x)], rel_ants)

        # get antenna positions for each baseline
        ant1_pos = np.array(map(lambda x: self.ant_pos[x[0]], str_bls))
        ant2_pos = np.array(map(lambda x: self.ant_pos[x[1]], str_bls))

        # get vector in u-v plane
        wavelengths = 2.9979e8 / (self.freqs * 1e6)
        self.uv_vecs = (ant2_pos - ant1_pos).reshape(Nbls, 3, -1) / wavelengths

        # Get direction unit vector, s-hat
        x, y, z = hp.pix2vec(self.beam_nside, np.arange(self.beam_npix))
        self.s_hat = np.array([y, x, z])

        # get phase map in topocentric frame (sqrt of full phase term)
        self.phase = np.exp( -1j * np.pi * np.einsum('ijk,jl->ikl', self.uv_vecs, self.s_hat) )

        # build beams for each antenna in topocentric coordinates
        self.build_beams(ant_inds=ant_inds, one_beam_model=one_beam_model)

        # create phase and beam product maps (Npol, Nbl, Nfreqs, Npix)
        self.beam_phs = self.ant_beam_models.reshape(self.Npols, -1, self.Nfreqs, self.beam_npix) * self.phase

        # loop over JDs
        self.vis_data = []
        for jd in JD_array:

            # map through each antenna and project polarized, multi-frequency beam model onto sky
            if one_beam_model == True:
                self.proj_beam_phs = self.project_beams(jd, self.sky_theta, self.sky_phi, beam_models=self.beam_phs, output=True, interp=interp)

            else:
                #proj_ant_beam_models = self.project_beams(JD, self.sky_theta, self.sky_phi,
                #beam_models=np.array(map(lambda x: self.ant_beam_models[str(x)], self.ant_nums)), output=True)
                self.proj_beam_phs = self.project_beams(jd, self.sky_theta, self.sky_phi, beam_models=self.beam_phs, output=True, interp=interp)

            # calculate visibility
            if pool is None:
                M = map
            else:
                M = pool.map

            if one_beam_model == True:
                if self.onepol is True:
                    vis = np.einsum("ijkl, ikl -> ijk", self.proj_beam_phs**2, self.sky_models)
                else:
                    vis = np.einsum("ijkl, imkl, mjkl -> imjk", self.proj_beam_phs, self.sky_models, self.proj_beam_phs)
            else:
                vis = np.array(map(lambda x: np.einsum("ijk, jk, ljk -> ilj", self.proj_beam_phs[rel_ant2ind[x[1][0]]],
                      self.sky_models, self.proj_beam_phs[rel_ant2ind[x[1][1]]]), enumerate(str_bls)))

            self.vis_data.append(vis)

        if self.onepol is True:
            self.vis_data = np.moveaxis(np.array(self.vis_data), 0, 2)
        else:
            self.vis_data = np.moveaxis(np.array(self.vis_data), 0, 3)
        
        # normalize by the effective-beam solid angle
        if one_beam_model == True:
            self.ant_beam_sa = np.repeat(np.einsum('ijk,ijk->ij', self.ant_beam_models, self.ant_beam_models)[ :, np.newaxis, :], Nbls, 1)
        else:
            self.ant_beam_sa = np.array(map(lambda x: np.einsum('ijk,ljk->ilj', self.proj_ant_beam_models[rel_ant2ind[x[1][0]]],
                                                            self.proj_ant_beam_models[rel_ant2ind[x[1][1]]]), enumerate(str_bls)))
        if self.onepol is True:
            self.vis_data /= self.ant_beam_sa[:, :, np.newaxis, :]
            self.vis_data_shape = "(1, Nbls, Ntimes, Nfreqs)"
        else:
            self.vis_data /= self.ant_beam_sa[:, :, :, np.newaxis, :]
            self.vis_data_shape = "(2, 2, Nbls, Ntimes, Nfreqs)"

        # write to file
        if write_miriad == True:

            # create uvdata object
            uvd = UVData()

            uvd.Nants_telescope    = self.Nants
            uvd.Nants_data         = len(np.unique(bl_array))
            uvd.Nbls               = Nbls
            uvd.Ntimes             = Ntimes
            uvd.Nblts              = Nbls * Ntimes
            uvd.Nfreqs             = self.Nfreqs
            uvd.Npols              = self.Nxpols
            uvd.Nspws              = 1
            uvd.antenna_numbers    = self.ant_nums
            uvd.antenna_names      = np.array(self.ant_nums, dtype=str)
            uvd.ant_1_array        = np.concatenate([map(lambda x: x[0], bl_array) for i in range(Ntimes)])
            uvd.ant_2_array        = np.concatenate([map(lambda x: x[1], bl_array) for i in range(Ntimes)])
            uvd.baseline_array     = np.concatenate([np.arange(Nbls) for i in range(Ntimes)])
            uvd.freq_array         = (self.freqs * 1e6).reshape(1, -1)
            uvd.time_array         = np.repeat(np.array(JD_array)[:, np.newaxis], Nbls, axis=1).ravel()
            uvd.channel_width      = (self.freqs * 1e6)[1] - (self.freqs * 1e6)[0]
            uvd.data_array         = self.vis_data.reshape(uvd.Nblts, 1, self.Nfreqs, self.Nxpols)
            uvd.flag_array         = np.ones_like(uvd.data_array, dtype=np.bool)
            uvd.history            = " "
            uvd.instrument         = " "
            uvd.integration_time   = 10.7
            uvd.lst_array          = np.repeat(np.array(map(lambda x: self.JD2LST(self.loc, x), JD_array))[:, np.newaxis], Nbls, axis=1).ravel()
            uvd.nsample_array      = np.ones_like(uvd.data_array, dtype=np.float)
            uvd.object_name        = "GSM"
            uvd.phase_type         = "drift"
            uvd.polarization_array = np.array(map(lambda x: {"XX":-5,"YY":-6,"XY":-7,"YX":-8}[x], self.xpols))
            uvd.spw_array          = np.array([0])
            uvd.telescope_location = np.array([ 6378137.*np.cos(self.loc.lon)*np.cos(self.loc.lat),
                                                6378137.*np.sin(self.loc.lon)*np.cos(self.loc.lat),
                                                6378137.*np.sin(self.loc.lat) ])
            uvd.telescope_name     = " "
            uvd.uvw_array          = np.ones((uvd.Nblts, 3))
            uvd.vis_units          = "Jy"
            zen_dec, zen_ra        = self.loc.radec_of(0, np.pi/2)
            uvd.zenith_dec         = np.ones(uvd.Nblts) * zen_dec
            uvd.zenith_ra          = np.ones(uvd.Nblts) * zen_ra

            uvd.write_miriad(fname, clobber=clobber)


































