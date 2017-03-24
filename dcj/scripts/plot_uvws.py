#!/usr/bin/env python
import aipy as a, numpy as np
import optparse, sys, os
from matplotlib.pyplot import *

o = optparse.OptionParser()
o.set_usage('plot_uvws.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--uvdata',action='store_true',
    help='Use uvdata')
o.add_option('--phase',action='store_true',
    help='Rotate to the ra-dec of the first pointing of the drift scan (only works in uvdata mode on miriad files).')
#a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

uvws = []
phased_uvws = []
if opts.uvdata:
    import uvdata
    for uvfile in args:
        UV = uvdata.uv.UVData()
        if uvfile.endswith('uvfits'):
            UV.read_uvfits(uvfile)
        elif uvfile.startswith('zen'):
            UV.read_miriad(uvfile)
        uvws.append(UV.uvw_array.copy())
        if opts.phase:
            UV.phase(time = UV.time_array[0])
            phased_uvws.append(UV.uvw_array)
    uvws = np.concatenate(uvws).T
    if opts.phase:
        phased_uvws = np.concatenate(phased_uvws).T
        
else:
    
    for uvfile in args:
        uv = a.miriad.UV(uvfile)
        cur_time = 0
        for (uvw,t,(i,j)),d in uv.all():
            uvws.append(uvw.copy())
            cur_time=t
    uvws = np.array(uvws)
    


uvw_lengths = np.array([np.sqrt(np.dot(uvw,uvw)) for uvw in uvws])

print "read {n} uvws".format(n=len(uvws))
print "uvw range (min, max):",np.min(uvw_lengths[uvw_lengths>0]),np.max(uvw_lengths)
sys.stdout.flush()
subplot(131)
plot(uvws[:,0],uvws[:,1],'o')
if opts.phase:
    plot(phased_uvws[:,0],phased_uvws[:,1],'.')
xlabel('u')
ylabel('v')
subplot(132)
plot(uvws[:,0],uvws[:,2],'o')
if opts.phase:
    plot(phased_uvws[:,0],phased_uvws[:,2],'.')
xlabel('u')
ylabel('w')
subplot(133)
plot(np.sort(uvw_lengths[uvw_lengths>0]))

show()
