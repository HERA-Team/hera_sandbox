#! /usr/bin/env python

#import matplotlib as mpl
#mpl.rcParams['font.size']  = 18
#mpl.rcParams['legend.numpoints']  = 1
#mpl.rcParams['legend.frameon'] = False
#mpl.rcParams['legend.fontsize'] = 8
#mpl.rcParams['figure.dpi'] = 400
#mpl.rcParams['savefig.dpi'] = 400
#mpl.rcParams['savefig.format'] ='png'
#mpl.rcParams['lines.markeredgewidth'] = 0
#mpl.rcParams['lines.markersize'] = 7
#mpl.rcParams['lines.linewidth'] = 2
import aipy as a, capo.arp as arp, capo.frf_conv as fringe, capo.zsa as zsa
import capo, numpy as n, pylab as p, sys, os, optparse
import matplotlib.cm as cmx
import matplotlib.colors as colors
from scipy.special import erf
import scipy.stats as stats

import pickle

def skew(cenwid, bins):
        return n.exp(-(bins-cenwid[0])**2/(2*cenwid[1]**2))*(1+erf(cenwid[2]*(bins-cenwid[0])/(n.sqrt(2)*cenwid[1]))) 
"""
plot the frf response and print out the effective integration time, same inputs as frf_filter.py
needs as input one file for baselines etc
"""

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('--frpad',type='string',default='1.0',help='make the fringe rate convolution longer by this factor, can accept comma separated list (default 1.0)')
o.add_option('--seps',type=str,
    help='list of seperations to use, ex 0,1;-1,1')
o.add_option('--plot',action='store_true',
        help='Plot the Fringe Rate Width')
o.add_option('--maxfr',type='float', default='1.2',
        help='max fringe rate for filter')
o.add_option('--inttime',type='float', default='50',
        help='set the inttime used to calculate the fringe-rate-filter')
o.add_option('--lstbintime',type=float,default=42.949,
        help='set the integration time of the actual data')
o.add_option('--output',action='store_true',
        help='output each frf as a pkl that can be plotted by plot_frf')
o.add_option('--noise',type='float',default=0.0,
        help='Supply noise filter level in jansky')
o.add_option('--mychan', type='float',default=0,
        help='manually supply channel to make frf, defualts to FLOOR(nchan/2)')
opts,args = o.parse_args(sys.argv[1:])

def get_colors(N):
    '''Returns function with N unique colors'''
    norm=colors.Normalize(vmin=0,vmax=N-1)
    scal_map=cmx.ScalarMappable(norm=norm,cmap='hsv')
    def map_index_to_rgb(index):
        return scal_map.to_rgba(index)
    return map_index_to_rgb

def noise(size):
    return n.random.random(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))
def noise_equivalent_bandwidth(t,H):
    #return 1/n.max(n.abs(H))**2 * n.sum(n.abs(H)**2)*n.diff(t)[0]
    #peak normalize
    H /= n.abs(H).max()
    return n.sum(n.abs(H)**2)*n.diff(t)[0]
    #return n.sum(H*H.conj())/(n.sum(H)*n.sum(H.conj()))


freqs = n.linspace(0.1,0.2,num=203)
aa = a.cal.get_aa(opts.cal, freqs)
nchan = len(freqs)
jy2T = capo.pspec.jy2T(freqs)
bls1,conj = capo.red.group_redundant_bls(aa.ant_layout)
NOISE= opts.noise
#pol = a.miriad.pol2str[uv['pol']]

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
if opts.mychan ==0:
    mychan = n.floor(nchan/2)
else: mychan = opts.mychan
#mychan = 160

#manually calculate the inttime so frf time bins match data
#(uvw,t1,(i,j)),d = uv.read()
#(uvw,t2,(i,j)),d = uv.read()
#while t1 == t2:
#    (uvw,t2,(i,j)),d = uv.read()
#inttime = (t2-t1)* (3600*24)
#del(uv)


#set channel to make frf
#mychan=101

##use calculated inttime to generate correct frf bins
frbins = fringe.gen_frbins(opts.inttime)
#frbins = n.arange( -.5/opts.inttime+5e-5/2, .5/opts.inttime,5e-5)
#DEFAULT_FRBINS = n.arange(-.01+5e-5/2,.01,5e-5) # Hz


print 'These are the separations that we are going to use ', seps
print "calculating fringe profile at channel ",mychan
#Get the fir filters for the separation used.
fig_firs,ax_firs=p.subplots(1)
fig_frp,ax_frp=p.subplots(1)
noise_array = {}
noise_array_top_hat = {}
noise_array_pre_filter = {}
for cnt,pad in enumerate(frpads):
    firs = {}
    frps = {}
    noise_array[pad] = {}
    noise_array_top_hat[pad] = {}
    noise_array_pre_filter[pad] = {}
    for sep in seps:
        c = 0 
        while c != -1:
            ij = map(int, sep2ij[sep].split(',')[c].split('_'))
            bl = a.miriad.ij2bl(*ij)
            if blconj[bl]: c+=1
            else: break
        frp, bins = fringe.aa_to_fr_profile(aa, ij, mychan,bins=frbins,frpad=pad, pol=opts.pol)
        timebins, firs[sep] = fringe.frp_to_firs(frp, bins, aa.get_afreqs(), fq0=aa.get_afreqs()[mychan], limit_xtalk=True,frpad=pad)


        if False and pad ==1:
            delta=prms0[-1]/n.sqrt(1+prms0[-1]**2)
            print 'model fit parameters: ',prms0
            print 'norm is: ', n.sum(frp)
            print 'mean is: ', n.sum(bins*frp)/n.sum(frp)
            mn= n.sum(bins*frp)/n.sum(frp)
            sq= n.sqrt(n.sum((bins-mn)**2*frp)/n.sum(frp))
            sk= n.sum(((bins-mn)/sq)**3*frp)/n.sum(frp)
            ftsk= (4-n.pi)/2.* (delta*n.sqrt(2/n.pi))**3/(1-2*delta**2/n.pi)**(1.5)
            print 'actual skew is: ', sk
            print 'fitted skew is: ', ftsk

        frps[sep], frp_freqs = fringe.fir_to_frp(firs[sep],
                tbins=timebins*opts.lstbintime/opts.inttime)
        baselines = ''.join(sep2ij[sep] for sep in seps)

    #print the effective integration time for each seperation
    for sep in seps:
        print "sep",sep
        print "   NEBW T_eff = ",
        print 1.2/noise_equivalent_bandwidth(frp_freqs,frps[sep][mychan])

        if NOISE > 0.: # Create a fake EoR signal to inject
            print 'FILTERING WHITE NOSIE at {0} Jy ...'.format(NOISE) ,
            sys.stdout.flush()
            ## this last term is to make the power spectrum equal 
            ## to expected noise line if 21cmSense noise only no filter is run
            ed = noise((nchan,711))  * NOISE 
            ed1 = n.copy(ed)
            noise_array_pre_filter[pad][sep] = n.copy(ed)
            top_hat = n.ones(n.int(n.round(dt_50))/opts.inttime) / ( n.int(n.round(dt_50))/opts.inttime)

            for cnt,ch in enumerate(n.arange(nchan)):
                ed[cnt] = n.convolve(ed[cnt], n.conj(firs[sep][ch]), mode='same')
                num = n.ones_like(ed[cnt].real)
                w = n.convolve ( num , n.abs(firs[sep][ch]), mode='same')
                ed[cnt] /= w
                ed1[cnt] = n.convolve(ed1[cnt], top_hat, mode='same')

            if blconj[bl]: 
                ed = n.conj(ed)
                ed1 = n.conj(ed1)
            noise_array[pad][sep] = n.transpose( ed.T *jy2T, [1,0])
            noise_array_top_hat[pad][sep] = n.transpose( ed1.T *jy2T, [1,0])

            if PLOT and True:
               fig,axes = p.subplots(nrows=3,ncols=1)
               p.subplot(311); p.plot( noise_array_pre_filter[pad][sep][mychan].real,alpha=.5); #p.colorbar();
               p.subplot(311); p.plot( noise_array[pad][sep][mychan].real,'k-'); #p.colorbar(); p.show()
               p.subplot(311); p.plot( noise_array_top_hat[pad][sep][mychan].real,'r-'); #p.colorbar(); p.show()
               p.ylabel('Real')
               p.subplot(312); p.plot( noise_array_pre_filter[pad][sep][mychan].imag,alpha=.5); #p.colorbar();
               p.subplot(312); p.plot( noise_array[pad][sep][mychan].imag,'k-'); #p.colorbar(); p.show()
               p.subplot(312); p.plot( noise_array_top_hat[pad][sep][mychan].imag,'r-'); #p.colorbar(); p.show()
               p.ylabel('Imag')
               p.subplot(313); p.plot( n.angle(noise_array_pre_filter[pad][sep][mychan]),alpha=.5); #p.colorbar();
               p.subplot(313); p.plot( n.angle(noise_array[pad][sep][mychan]),'k-'); #p.colorbar(); p.show()
               p.subplot(313); p.plot( n.angle(noise_array_top_hat[pad][sep][mychan]),'r-'); #p.colorbar(); p.show()
               p.ylabel('Phase')
               p.suptitle('%d_%d'%a.miriad.bl2ij(bl))
            print '[Done]'
    if opts.output: 
        filename = 'fringe_rate_profile_pad{pad}_int{inttime}_chan{chan:d}_pol{pol}.pkl'.format(pad=pad,inttime=opts.inttime,chan=int(mychan),pol=opts.pol)
        print("saving to: {filename}".format(filename=filename))
        f = open(filename,'w')
        frps['chan'] = mychan
        frps['freqs'] = frp_freqs
        frps['frp'] = frp
        frps['firs'] = firs
        frps['timebins'] = timebins
        pickle.dump(frps,f)
        f.close()
if PLOT:
    ax_frp.plot(frp_freqs*1e3,frp,label='frp0',color='black')
    #ax_frp.plot(frp_freqs*1e3,skew([mn,sq,sk],bins),'k--',label='data')
    ax_frp.set_title('Fitted Fringe Rate Profile')
    ax_frp.set_xlabel('Fringe Rate [mili Hz]')
    ax_frp.set_xlim([-.7,1.5])
    ax_frp.set_ylim([0,1])
    ax_frp.legend(loc='best')

    ax_firs.set_title('Fringe Rates Filter Widths')
    ax_firs.legend(loc='best')
    p.show()

#f= open('frf_diagnose_parsons_42.pkl','w')


#pickle.dump(frps,f)

#f.close()
