#!/usr/bin/python
"""
UN-Scale data by MULTIPLYING by spectrum nu^alpha
"""
import aipy, numpy as np, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('freqscale.py [options] *.uv')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True)
o.add_option('--alpha',dest='alpha',type='float',default=-1.5,help='Frequency scaling')
opts,args = o.parse_args(sys.argv[1:])

uv = aipy.miriad.UV(args[0])
aa = aipy.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

freq = aa.get_freqs()
sc=np.power(freq/0.15,opts.alpha)
        
for filename in args:
    outfile = filename+'u'
    print filename,'->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping'
        continue
    
    def mfunc(uv,p,d,f):
        d = d*sc
        return p,d,f
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,raw=True,append2hist='\n FREQSCALE: nu^%f \n'%opts.alpha)
    del uvo,uvi
