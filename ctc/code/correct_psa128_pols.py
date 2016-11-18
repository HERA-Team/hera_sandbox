#!/usr/bin/env python

import aipy
import capo
import numpy
import sys, optparse, os

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, ant=True)
o.add_option('-a','--ant',default=[], help='Comma-separated antenna numbers to swap pols.')
opts,args = o.parse_args(sys.argv[1:])

ANTS = []
for a in opts.ant.split(','):
    ANTS.append(int(a)) #antennas that need to be swapped x->y, y->x

def mfunc(uv,p,d,f):
    ant1,ant2 = p[2]
    if ant1 not in ANTS and ant2 not in ANTS: return p,d,f
    if ant1 in ANTS:
        if pol[0] == 'x': newpol1 = 'y'
        if pol[0] == 'y': newpol1 = 'x'
    else: newpol1 = pol[0]
    if ant2 in ANTS:
        if pol[1] == 'x': newpol2 = 'y'
        if pol[1] == 'y': newpol2 = 'x'
    else: newpol2 = pol[1]
    newpol = newpol1+newpol2
    index = numpy.where(t[newpol]['times'] == p[1])[0][0] #XXX times must match exactly
    d = data[newpol][p[2]][newpol][index] #collect information from correct pol file
    f = flags[newpol][p[2]][newpol][index]
    return p,d,f

pols = ['xx','xy','yx','yy']
for filename in args:
    files,t,data,flags = {},{},{},{}
    for pol in pols:
        files[pol] = filename.replace('xx',pol) #dictionary of files by pol
        t[pol],data[pol],flags[pol] = capo.miriad.read_files([files[pol]],antstr='all',polstr=pol,verbose=True) #read all the pol files
    for pol in files: #loop through 4 pols
        uvi = aipy.miriad.UV(files[pol])
        print files[pol], '->', files[pol]+'c'
        if os.path.exists(files[pol]+'c'):
            print '   File exists... skipping.'
            continue
        uvo = aipy.miriad.UV(files[pol]+'c',status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi,raw=True,mfunc=mfunc,append2hist='CORRECT_PSA128_POLS:'+' '.join(sys.argv)+'\n')
        
