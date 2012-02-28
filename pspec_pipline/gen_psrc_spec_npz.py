#! /usr/bin/env python
import numpy as n, aipy as a
import fg_sim

rfi_file = n.load('rfi_flags.npz')
freqs = rfi_file['freqs']
aa = a.cal.get_aa('psa746_v008', freqs)
POL = 'xx'

umag150 = 32.
         
umags = n.array([16., 32., 64., 128.])
specs = []
for umag150 in umags:
    print umag150
    bl_ns = umag150 / .150
    #specs.append(fg_sim.sim_srcs(aa, pol=POL, bl_ns=bl_ns, lo_cutoff_Jy=.1, hi_cutoff_Jy=5., std_index=.25))
    specs.append(fg_sim.sim_srcs(aa, pol=POL, bl_ns=bl_ns, lo_cutoff_Jy=.1, hi_cutoff_Jy=100., std_index=.25))

dat = {'freqs': freqs, 'umags': umags}
dat['spec'] = n.array(specs)

filename = 'psrc_spec.npz'
print 'Writing', filename
n.savez(filename, **dat)
