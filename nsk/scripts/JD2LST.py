#!/usr/bin/env python2.7

from astropy.time import Time
import argparse

a = argparse.ArgumentParser()
a.add_argument('--jd', default=[], type=float, help='JD float to convert to LST')
a.add_argument('--lon', default=[], type=float, help='longitude float of observer in degrees East')

def JD2LST(JD, longitude):
    """
    Input:
    JD : float, julian date
    longitude : float, longitude in degrees East

    Output:
    Local Apparent Sidreal Time in Hour Angle
    """
    if type(JD):
        pass
    t = Time(JD, format='jd')
    return t.sidereal_time('apparent', longitude=longitude).value


if __name__ == "__main__":
    args = a.parse_args()
    if args.jd is not None and args.lon is not None:
        print JD2LST(args.jd, args.lon)
