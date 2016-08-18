#! /usr/bin/env python
from astropy import time
import sys

print time.Time(float(sys.argv[-1]),scale='utc',format='gps').iso
