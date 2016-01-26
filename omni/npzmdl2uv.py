#! /usr/bin/env python

import aipy
import numpy
import capo
import os,sys
import optparse

##### 
# Outputs a UV file containing the model visibilities outputted by Omnical
#####


### Options ###
o = optparse.OptionParser()
o.set_usage('npzmdl2uv.py *npz')
aipy.scripting.add_standard_options(o,pol=True)
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

#Set up Miriad file
def mkuv(filename,data):
    uv = aipy.miriad.UV(filename, status='new')
    uv._wrhd('obstype','croscorrelation')
    uv._wrhd('history','Omnical Model: created file.\nOmnical Model: ' + ' '.join(sys.argv) + '\n')
    uv.add_var('telescop' ,'a'); uv['telescop'] = 'PAPER128'
    uv.add_var('operator' ,'a'); uv['operator'] = 'PAPER128'
    uv.add_var('version' ,'a'); uv['version'] = '1.3.0'
    uv.add_var('epoch' ,'r'); uv['epoch'] = 2000.
    uv.add_var('source'  ,'a'); uv['source'] = 'zenith'
    uv.add_var('nchan' ,'i'); uv['nchan'] = 203
    uv.add_var('sdf' ,'d'); uv['sdf'] = (0.2-0.1)/203
    uv.add_var('sfreq' ,'d'); uv['sfreq'] = 0.1
    uv.add_var('freq' ,'d'); uv['freq'] = 0.1
    uv.add_var('restfreq' ,'d'); uv['restfreq'] = 0.1
    uv.add_var('nschan' ,'i'); uv['nschan'] = uv['nchan']
    uv.add_var('npol' ,'i'); uv['npol'] = 1
    uv.add_var('nspect' ,'i'); uv['nspect'] = 1
    uv.add_var('nants' ,'i'); uv['nants'] = 128
    uv.add_var('vsource' ,'r'); uv['vsource'] = 0.
    uv.add_var('ischan'  ,'i'); uv['ischan'] = 1
    uv.add_var('tscale'  ,'r'); uv['tscale'] = 0.
    uv.add_var('veldop'  ,'r'); uv['veldop'] = 0.
    uv.add_var('coord' ,'d')
    uv.add_var('time' ,'d')
    uv.add_var('lst' ,'d')
    uv.add_var('ra' ,'d')
    uv.add_var('obsra' ,'d')
    uv.add_var('baseline' ,'r')
    uv.add_var('pol' ,'i')
    uv.add_var('inttime' ,'r')
    uv.add_var('latitud' ,'d'); uv['latitud'] = aa.lat
    uv.add_var('dec' ,'d'); uv['dec'] = aa.lat
    uv.add_var('obsdec' ,'d'); uv['obsdec'] = aa.lat
    uv.add_var('longitu' ,'d'); uv['longitu'] = aa.long
    uv.add_var('antpos' ,'d'); uv['antpos'] = numpy.array([ant.pos for ant in aa], dtype=numpy.double).transpose().flatten() #transpose is miriad convention

    for bl in data.keys():
        uv['inttime'] = ((max(data[bl].keys())-min(data[bl].keys()))/(len(data[bl].keys())-1))*24*3600
        times = numpy.sort(data[bl].keys())
        for t,time in enumerate(times):
            i,j = int(bl.split('_')[0]),int(bl.split('_')[1])
            bsln = aa.get_baseline(i,j)
            preamble = (bsln, time, (i,j))  
            aa.set_jultime(time)
            uv['time'] = time
            uv['lst'] = aa.sidereal_time()
            uv['ra'] = aa.sidereal_time()
            uv['obsra'] = aa.sidereal_time()
            uv['pol'] = aipy.miriad.str2pol[opts.pol]
            uv['coord'] = bsln
            uv['baseline'] = float(aipy.miriad.ij2bl(i,j))
            flags = numpy.ma.masked_where(data[bl][time]==0,data[bl][time]).mask
            finaldata = numpy.ma.masked_array(data[bl][time],flags)
            uv.write(preamble, finaldata, flags)    
    del(uv)


#Read NPZ files
aa = aipy.cal.get_aa('psa6622_v003', 0.001,0.1,203)
for file in args:
    data = {}
    npz = numpy.load(file)
    jds = npz['jds']
    for key in npz.keys():
        if key[0] == '<' and key[-2:] == opts.pol:
            vis = npz[key] 
            a1 = key.split(' ')[0].split(',')[0].split('<')[1]
            a2 = key.split(' ')[0].split(',')[1].split('>')[0]
            bl = a1+'_'+a2
            data[bl] = {}
            for t in range(len(vis)):
                data[bl][jds[t]] = vis[t] #data saved by bl and time
    name = file.split('npz')[0]+opts.pol+'.uvcRRE'
    print file, '->', name+'M'
    if os.path.exists(name+'M'):
        print '    File exists... skipping.'
        continue
    mkuv(name+'M',data)

 
