#! /usr/bin/env python
import healpy as hp
import numpy as n
import optparse,sys,os
"""
smooth_beam.py --fwhm=10 <files>
smooth healpix beam files (also works on anything else that makes sense in log space (no negatives)
"""

o = optparse.OptionParser()
o.add_option('--fwhm',default=10.0,type='float',
    help='output resolution [deg]')
o.add_option('--clip',default=0,type='float',
    help="""wacky outliers can bogus a smooth. Throw out things above this many sigma (does the operation on log data)
    [default None]""")
opts,args = o.parse_args(sys.argv[1:])
for file in args:
    outfile = file[:-5]+'_smooth.fits'
    print file,'->',outfile
    if os.path.exists(outfile):
        print "file exists"
        continue
    beam = hp.read_map(file)
    beam_rms = n.std(n.log10(beam[beam!=0]))
    beam = n.ma.masked_where(n.abs(beam)>10**(opts.clip*beam_rms),beam)
    beam_smooth = hp.smoothing(beam.filled(0),fwhm=opts.fwhm,degree=True)
    hp.write_map(outfile,beam_smooth)
