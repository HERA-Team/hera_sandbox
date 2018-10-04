from __future__ import absolute_import, division, print_function

import numpy as np
import collections
import six
import warnings
from astropy.time import Time
from astropy import coordinates as crd
from astropy import units as unt
import os


def RA2Time(ra, anchor_jd, latitude=-30.72, longitude=21.42, return_lst=True, tolerance=1e-5, maxiter=5):
    """
    Given a right ascension coordinate in J2000 (ICRS), return the
    Local Sidereal Time (LST) or Julian Date on the day of anchor_jd given
    an observer at the specified latitude and longitude on Earth.

    Parameters
    ----------
    ra : float
        Right Ascension in J2000 (ICRS) frame.

    anchor_jd : float
        Julian Date to anchor the conversion from ra to time, because RA wraps with
        respect to time.

    latitude : float
        Latitude on Earth of observer in degrees.

    longitude : float
        Longitude on Earth of observer in degrees.

    return_lst : bool
        If True, return time as local sidereal time in [radian], else return as julian date.

    tolerance : float
        Convergence tolerance in Julian Date for minimizer.

    maxiter : int
        Maximum number of iterations in minimization.
    
    Returns
    -------
    time : float
        LST [radian] or JD when ra is closest to zenith of an observer
    """
    # setup Earth Location
    loc = crd.EarthLocation(lat=latitude * unt.deg, lon=longitude * unt.deg)

    # Enter loop
    k = 0
    time = anchor_jd
    while True:
        # setup Time object with initial Julian Date and a step forward
        t1 = Time(time, format='jd', scale='utc')
        t2 = Time(time+1e-4, format='jd', scale='utc')

        # calculate derivative of d_ra / d_jd
        ra1 = crd.SkyCoord(frame='altaz', alt=90 * unt.deg, az=0 * unt.deg, location=loc, obstime=t1).icrs.ra.degree
        ra2 = crd.SkyCoord(frame='altaz', alt=90 * unt.deg, az=0 * unt.deg, location=loc, obstime=t2).icrs.ra.degree
        d_ra_d_jd = (ra2 - ra1) / 1e-4

        # calculate d_jd from initial time to get desired RA
        if ra < ra1 - 180:
            ra1 -= 360
        d_jd = (ra - ra1) / d_ra_d_jd

        # get new time
        time = time + d_jd

        if np.abs(d_jd) <= tolerance:
            break
        elif k >= maxiter:
            break

        k += 1

    # convert to LST if desired
    if return_lst:
        t = Time(time, format='jd', scale='utc')
        time = t.sidereal_time('apparent', longitude=longitude * unt.deg).radian

    return time

