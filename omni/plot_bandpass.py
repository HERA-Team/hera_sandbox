#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys
import matplotlib.pyplot as plt

o = optparse.OptionParser()
o.set_usage('plot_bandpass.py [options] *uvcRRE.ucalbandpass.npz')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,pol=True)
opts,args = o.parse_args(sys.argv[1:])

for file in args:
    print 'Reading', file
    npz = numpy.load(file)
    bp = npz['bandpass']
    plt.plot(bp,label=file.split('.')[1])
plt.legend(prop={'size':8})
plt.show()    
