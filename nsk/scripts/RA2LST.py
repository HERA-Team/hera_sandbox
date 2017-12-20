"""
given an RA and a longitude on Earth,
calculate LST when that RA is at zenith.
"""
import ephem
import numpy as np
import datetime
import pytz

def RA2LST(RA, lon):
    """
    RA : float
        right ascension (J2000) in degrees

    lon : float
        longitude East of observer in degrees

    return LST_RA (LST in hours)
    """
    # get observer
    obs = ephem.Observer()
    obs.lon = lon * np.pi / 180.0
    obs.date = datetime.datetime(2000, 03, 20, 12, 0, 0, 0, pytz.UTC)

    # get RA at zenith of observer in degrees
    ra_now = obs.radec_of(0, np.pi/2)[0] * 180 / np.pi

    # get LST of observer
    LST_now = obs.sidereal_time() * 12.0 / np.pi 

    # get the LST of the RA via difference
    LST_RA = LST_now + (RA - ra_now) / 15.0
    return LST_RA





