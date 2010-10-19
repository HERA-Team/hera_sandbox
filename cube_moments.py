#!/usr/bin/env python
#
#  cube_moments.py
#  
#
#  Created by Danny Jacobs on 1/15/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
"""
Compress input fits cube into moment maps
cube_moments.py *.fits 
for each input <file>.fits computes
<file>.m0.fits = average(flux) 
<file>.m1.fits = average(diff(log(flux)))  ~= spectral ind
<file>.var.fits = sqrt(var(flux)) = sqrt(average(flux^2)-average(flux)^2)
    All along frequency dimension

"""
o = optparse.OptionParser()
a.scripting.add_standard_options(o)
opts, args = o.parse_args(sys.argv[1:])

def append_file(filename,append):
    return '.'.join(file.split('.')[:-1])+append+'.fits'

for file in args:
    print file
    data,kwds = a.img.from_fits(file)
    #mask off the NaN
    data = n.ma.array(data,mask=n.isnan(data))
    for i,v in enumerate(kwds['axes']):
        if v.lower().startswith('freq'):
            ax=i
    
    m0 = n.array(n.ma.average(data,axis=ax))
    m0.shape = m0.shape + (1,)
    
    m1 = n.array(n.ma.average(n.diff(n.log10(data),axis=ax),axis=ax))
    m1.shape = m1.shape + (1,)
    var = n.array(n.ma.sqrt(n.var(data,axis=ax)))
    var.shape = var.shape + (1,)
    print m0.shape,n.max(m0),n.min(m0)    
    #write out files
    outfile = append_file(file,'m0')
    print '\t > '+ outfile
    a.img.from_fits_to_fits(file,outfile,m0,kwds,
         ' '.join(args) + "\n" + "computed zeroth moment in freq direction")

    outfile = append_file(file,'m1')
    print '\t > '+ outfile
    a.img.from_fits_to_fits(file,outfile,m1,kwds,
        ' '.join(args) + "\n" + "computed first moment in freq direction")
    
    outfile = append_file(file,'var')
    print '\t > '+ outfile
    a.img.from_fits_to_fits(file,outfile,var,kwds,
        ' '.join(args) + "\n" + "computed rms in freq direction")