"""
beam_sim.py
----------

Routines for handling beam files in healpix

Nick Kern
Sept. 2017
"""

# Load Modules
import matplotlib
matplotlib.use('Agg')
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

    def JD2LST(self, loc, JD):
        """
        use ephem to convert from JD to LST
        """
        loc2 = loc.copy()
        loc2.date = Time(JD, format="jd").datetime
        return loc2.sidereal_time()

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

class Beam_Model(Helper):
    """
    Beam_Model
    ----------

    class for handling of beam models
    """
    def __init__(self, beamfile, loc=None, pols=None, freqs=None, pool=None, verbose=False, mask=True):
        """
        Load and configure beam models

        Input:
        ------
        beamfile : str
            pyuvdata beamfits file

        pols : list, default=['X']
            list of polarizations of antenna feeds
            to put in self.beam_models
            options=['X', 'Y'], for North and/or East
            North (Y) is assumed to be 90 deg rotation of East (X)

        freqs : ndarray, default=None
            frequency channels in MHz, defaults to
            channels in steps of 1 MHz

        beam_data_path : str
        	path to beam data in fits format with 2 HDUs
        	with healpix beam models in 0th hdu w/ shape=(Npix, Nfreq)
        	and beam frequencies in 1st hdu w/ shape=(Nfreq,)

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

        # select polarization
        if pols is not None:
            if len(pol_arr) == 1:
                self.print_message("...only polarization in beamfits is %s"%pol_arr, verbose=verbose)
            else:
                try:
                    pol_select  = map(lambda x: {"X":0, "Y":1}[x], pols)
                    beam_models = beam_models[pol_select]
                except KeyError:
                    raise Exception("...didn't recognize pol %s, try 'X' or 'Y'"%pols)
                self.print_message("...using pol(s) %s"%pols, verbose=verbose)

        # 1D interpolation of frequency axis if desired
        if freqs is not None:
            beam_models = interpolate.interp1d(beam_freqs, beam_models, axis=1,\
                                    bounds_error=False, fill_value='extrapolate')(freqs)
            beam_freqs = freqs

        # Get theta and phi arrays
        self.beam_nside = uvb.nside
        self.beam_npix = uvb.Npixels
        beam_theta, beam_phi = hp.pix2ang(self.beam_nside, np.arange(self.beam_npix), lonlat=False)

        # mask beam models below horizon
        mask = (beam_phi > np.pi/2.0) & (beam_phi < 3*np.pi/2.0)
        beam_models[:, :, mask] - 0.0

        # assign vars to class
        self.beam_models = beam_models
        self.beam_freqs  = beam_freqs
        self.beam_theta  = beam_theta
        self.beam_phi    = beam_phi

    def rotate_beam(self, beam_models, beam_theta, beam_phi, rot=[0,0,0]):
    	"""
    	rotate beam models

    	Input:
    	-------
    	beam_models : ndarray, shape=(Nfreq, Npix)
    		multi-frequency beam models in healpix format

    	beam_theta : ndarray, shape=(Npix,)
			healpix colatitude in radians

		beam_phi : ndarray, shape=(Npix,)
			healpix longitude in radians

		rot : array, default=[0,0,0], shape=(3,)
			rotation angles in radians (lon, lat, omega)

		pool : pool object, default=None
			multiprocessing pooling object
		"""
    	R = hp.Rotator(rot=rot, deg=False)
    	beam_theta2, beam_phi2 = R(beam_theta, beam_phi)
        beam_models = hp.ma(self.healpix_interp(beam_models, self.beam_nside, beam_theta2, beam_phi2))
    	return beam_models

    def project_beams(self, JD, sky_theta, sky_phi, beam_models=None, obs_lat=None, obs_lon=None,
    					freqs=None, output=False, pool=None):
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

        # Apply rotation to equatorial coordinates of observer
        hrot = hp.Rotator(rot=[obs_ra, obs_dec], coord=['G', 'C'], inv=True, deg=False)
        theta_eq, phi_eq = hrot(sky_theta, sky_phi)

        ## Interpolate beam at healpix values to get beam models
        ## projected onto the sky
        if pool is not None:
            M = pool.map
        else:
            M = map

        # Get direction unit vector, s-hat, the way it outputs, X==y, Y==z, and Z==x
        x, y, z = hp.pix2vec(self.beam_nside, np.arange(self.beam_npix))
        self.s_hat = np.array([y,z,x])
        self.s_hat = self.healpix_interp(self.s_hat, self.beam_nside, theta_eq, phi_eq).T

	    # iterate through polarization, then map through each beam model in frequency
        self.sky_beam_models = self.healpix_interp(beam_models, self.beam_nside, theta_eq, phi_eq)
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
    				cbar_kwargs={}):
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
