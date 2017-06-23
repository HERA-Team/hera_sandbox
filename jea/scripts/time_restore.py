import aipy as a
import sys
import capo
from glob import glob
import pylab as plt
import numpy as np
import time

import mytools as my
file = '/Users/jaguirre/Documents/PAPER/GlobalSignal/BaselinePulls/bl_0_46.npz'

t0 = my.stime('Reading npz of waterfall')
data = np.load(file)
my.etime(t0)

