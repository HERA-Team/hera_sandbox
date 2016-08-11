#! /usr/bin/env python
import numpy as n, aipy as a
from capo import omni
from capo import arp
import optparse,sys


o=optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv[1:])


