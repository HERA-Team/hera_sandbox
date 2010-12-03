#!/usr/bin/env /usr/lib64/casapy/bin/python
"""
CASA script for batch importing miriad files

D.Jacobs 2010
"""
import sys, os, optparse
#Setup casa environment (AOC)
sys.path.append('/usr/lib64/casapy/test/lib/python2.6')
sys.path.append('/usr/lib64/casapy/test/lib/python2.6/site-packages')
#setup miriad environment (AOC)
sys.path.append('/usr/local/miriad/linux64/bin/')
os.environ['MIRCAT']='/usr/local/miriad/cat/'
import casa

o = optparse.OptionParser()
o.set_usage('miriad2casams.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])




for file in args:
    fitsfile = file+'.fits'
    visfile = file+'.ms'
    print(file+' > '+fitsfile)
    os.system('fits in='+file+' op=uvout out='+fitsfile+' stokes=xx')
    print(fitsfile+' > ' + visfile)
    casa.importuvfits(fitsfile=fitsfile,vis=visfile)
