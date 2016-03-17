import numpy as np, capo, optparse, os, sys, aipy

"""
Quick hack of omni_check to return dictionary 
of gains for manipulation in IPython
"""

o = optparse.OptionParser()
o.set_usage('omni_get.py *.npz')
opts,args = o.parse_args(sys.argv[1:])

gains = {} 
for i, file in enumerate(args): #loop over files
    print 'Reading',file
    file = np.load(file)
    for key in file.keys(): #loop over antennas
        if key[0] != '<' and key[0] != '(' and key[0].isalpha() != True:
            gain = file[key]
            antnum = key[:-1]
            try: gains[antnum].append(gain)
            except KeyError: gains[antnum] = [gain]
            
for key in gains.keys():
    gains[key] = np.vstack(gains[key]) #cool thing to stack 2D arrays that only match in 1 dimension
    mk = np.ma.masked_where(gains[key] == 1,gains[key]).mask #flags
    gains[key] = np.ma.masked_array(gains[key],mask=mk) #masked array

from matplotlib import pylab
import IPython; IPython.embed()
