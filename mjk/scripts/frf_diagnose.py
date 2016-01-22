#! /usr/bin/env python
import aipy as a
import capo.arp as arp
import capo.frf_conv as fringe
import capo.zsa as zsa
import numpy as n, pylab as p
import sys, os, optparse
import matplotlib.cm as cmx
import matplotlib.colors as colors
"""
plot the frf response and print out the effective integration time, same inputs as frf_filter.py
needs as input one file for baselines etc
"""

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--frpad',type='string',default='1.0',help='make the fringe rate convolution longer by this factor, can accept comma separated list (default 1.0)')
o.add_option('--seps',type=str,
    help='list of seperations to use, ex 0,1;-1,1')
o.add_option('--plot',action='store_true',
        help='Plot the Fringe Rate Width')
opts,args = o.parse_args(sys.argv[1:])

def get_colors(N):
    '''Returns function with N unique colors'''
    norm=colors.Normalize(vmin=0,vmax=N-1)
    scal_map=cmx.ScalarMappable(norm=norm,cmap='hsv')
    def map_index_to_rgb(index):
        return scal_map.to_rgba(index)
    return map_index_to_rgb


uv = a.miriad.UV(args[0])
#nants = uv['nants']
#inttime = uv['inttime'] * 4 #integration time hack
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
nchan = uv['nchan']
#pol = a.miriad.pol2str[uv['pol']]
del(uv)

#Get only the antennas of interest
sep2ij, blconj, bl2sep = zsa.grid2ij(aa.ant_layout)

#print "Looking for baselines matching ", opts.ant
#ants = [ b[0] for b in a.scripting.parse_ants(opts.ant, nants) ]
#seps = [ bl2sep[b] for b in ants ]
#seps = n.unique(seps)
PLOT=opts.plot
seps = opts.seps.split(';')
frpads = [ float(x) for x in opts.frpad.split(',')]
npads=len(frpads)
cmap= get_colors(int(1.5*npads))
mychan = n.floor(nchan/2)
print 'These are the separations that we are going to use ', seps
print "calculating fringe profile at channel ",mychan
#Get the fir filters for the separation used.
for cnt,pad in enumerate(frpads):
    firs = {}
    for sep in seps:
        c = 0 
        while c != -1:
            ij = map(int, sep2ij[sep].split(',')[c].split('_'))
            bl = a.miriad.ij2bl(*ij)
            if blconj[bl]: c+=1
            else: break
        frp, bins = fringe.aa_to_fr_profile(aa, ij, mychan,frpad=pad)
        timebins, firs[sep] = fringe.frp_to_firs(frp, bins, aa.get_afreqs(), fq0=aa.get_afreqs()[mychan])
        
    baselines = ''.join(sep2ij[sep] for sep in seps)
    #times, data, flags = arp.get_dict_of_uv_data(args, baselines, pol, verbose=True)
    #lsts = [ aa.sidereal_time() for k in map(aa.set_jultime(), times) ]
    for sep in seps:
        if PLOT:
           # p.plot(timebins,firs[sep][mychan])
            p.plot(timebins,n.abs(firs[sep][mychan]),label='{0}'.format(pad),color=cmap(cnt))
            p.xlabel('s')
        envelope = n.abs(firs[sep][mychan])
        envelope /= n.max(envelope)
        dt = n.sqrt(n.sum(envelope*timebins**2)/n.sum(envelope))
        dt_50 = (timebins[envelope>0.5].max() - timebins[envelope>0.5].min())
        print "pad:", pad, "variance width ",sep, " [s]:",int(n.round(dt)),"50% width",int(n.round(dt_50))
if PLOT:
    p.title('Fringe Rates Filter widths')
    p.legend(loc='best')
    p.show()
