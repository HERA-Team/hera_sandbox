import capo.omni as omni
import numpy as n, pylab as p, capo.plot as plot
import sys


args = sys.argv()[1:]

m,g,v,x = omni.from_npz(args, verbose=True)
    
