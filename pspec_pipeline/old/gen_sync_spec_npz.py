#! /usr/bin/env python
import aipy as a, numpy as n
import fg_sim

rfi_file = n.load('rfi_flags.npz')
freqs = rfi_file['freqs']
hmap = a.map.Map(fromfits='haslam408.fits')
aa = a.cal.get_aa('psa746_v008', freqs)
POL = 'xx'
#JD = [2455746.8, 2455746.9, 2455747.0]
JD = 2455746.8
         
umags = n.array([16., 32, 64, 128])
specs = []
for umag150 in umags:
    print umag150
    bl_ns = umag150 / .150
    u,v = (bl_ns,0.)
    specs.append(fg_sim.sim_sync(freqs, hmap, aa, bl_ns=(u,v), pol=POL, jd=JD))

dat = {'freqs': freqs, 'umags': umags}
dat['spec'] = n.array(specs)

filename = 'sync_spec.npz'
print 'Writing', filename
n.savez(filename, **dat)
