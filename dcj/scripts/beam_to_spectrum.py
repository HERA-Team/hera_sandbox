#!/usr/bin/env python
#
#  beam_to_spectrum.py
#  
#
#  Created by Danny Jacobs on 3/4/10.
#  PAPER Project
#
"""
Average down a beamformed file to get a spectrum, as well as a catalog point.
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,glob,re

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True,pol=True)
opts, args = o.parse_args(sys.argv[1:])

def inc_filename(match):
    "match=ls explicit search criteria, eg name_*.txt, where * is replaced with next highest number startinga at 1" 
    files = glob.glob(match)
    if len(files)<1:
        return re.sub('\*','1',match)
    files = n.sort(files)
    i = re.search(re.escape(match.split('*')[0])+'(\d+)'+match.split('*')[1]
        ,files[-1]).groups(1)[0]
    return re.sub('\*',str(int(i)+1),match)
dbuf = []
dt = None
for file in args:
    outfile = inc_filename(opts.src+'_spectrum_*.txt')
    print file, '->', outfile
    uv = a.miriad.UV(file)
    freqs = n.arange(uv['sfreq'], uv['sfreq']+uv['nchan']*uv['sdf'], uv['sdf']) 
    a.scripting.uv_selector(uv, 'all', opts.pol)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if dt is None: dt = n.dtype([
                    ('time',n.float32),
                    ('spectrum',n.float32,d.shape),
                    ('flags',n.bool,f.shape)])
        dbuf.append((t,d,f))
dbuf = n.array(dbuf,dtype = dt)
_outfile = open(outfile,'w')
for line in uv['history'].split('\n'):
    _outfile.write('#%s\n'%(line))
_outfile.write('\n')

spec = n.ma.array(dbuf['spectrum'],mask=dbuf['flags'])
spectrum = n.vstack([
            freqs,
            n.ma.average(spec,axis=0),
            n.ma.sqrt(n.ma.var(spec,axis=0))]).transpose()
_outfile.write("#! name=%s \t nu=%3.5f \t S_nu=%3.2f \t S_nu_e=%3.2f\n"%(
        opts.src,
        n.ma.average(spectrum[:,0]),
        n.ma.average(spectrum[:,1]),
        n.ma.sqrt(n.ma.var(spectrum[:,2]))
        ))
_outfile.write('# nu [Ghz] \t S_nu [Jy] \t S_nu_e [Jy]\n')
for l in spectrum:
    _outfile.write("%3.6f \t %3.2f \t %3.2f\n" %( l[0],l[1],l[2]))
_outfile.close()

    