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
import aipy as a
import capo.arp as arp
import capo.frf_conv as fringe
import capo.zsa as zsa
import numpy as n, pylab as p
import sys, os, optparse
import matplotlib.cm as cmx
import matplotlib.colors as colors
from scipy.special import erf
import scipy.stats as stats

import pickle

def skew(cenwid, bins):
        return n.exp(-(bins-cenwid[0])**2/(2*cenwid[1]**2))*(1+erf(cenwid[2]*(bins-cenwid[0])/(n.sqrt(2)*cenwid[1]))) 


def noise(size):
    size = list(size)
    size[-1] *= 10
    size = tuple(size)
    return n.random.normal(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))

def gen_flags(size,f):
    size = list(size)
    size[-1] *= 10
    f_size = list(f.shape)
    start= int( ( size[-1] - f_size[-1] )/2.)
    end = int( ( size[-1] + f_size[-1])/2. )
    diff = end-start - f_size[-1]
    if diff != 0: end += diff
    ind = n.arange(start, end,1)
    flags = n.zeros(size)
    for ch in xrange(size[0]): flags[ch,ind] = n.copy(f[ch])
    return flags

def clip_array(d,size,axis=0):
    if size is None:
        print 'Must give valid size for clipping'
        print 'Returning origian larray'
        return d
    d_size = d.shape[axis]
    if d_size < size:
        print 'Array Must have large size than desired clip'
        print 'Clipping failed, returning original array'
        return d
    start= int( (d_size - size)/2.)
    end = int( (d_size + size)/2. )
    diff = size - ( end - start )
    if diff != 0: end += diff
    d = n.swapaxes(d,axis,0)
    _d = d[start:end]
    _d = n.swapaxes(_d,0, axis)
    return _d
"""
plot the frf response and print out the effective integration time, same inputs as frf_filter.py
needs as input one file for baselines etc
"""

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,pol=True)
o.add_option('--bl_scale',type='string',default='1.0',help='makes the baseline used in fringe rate filter creation longer by this factor, can accept comma separated list (default 1.0)')
o.add_option('--fr_width_scale', type='string', default='1.0',
        help='Artificially inflates width of Fringe Rate Filter by scale factor, can accept comma separated list (default 1.0)')
o.add_option('--seps',type=str,
    help='list of seperations to use, ex 0,1;-1,1')
o.add_option('--plot',action='store_true',
        help='Plot the Fringe Rate Width')
o.add_option('--chan',type='int',
        help='Channel number to calculate FRP')
o.add_option('--maxfr',type='float', default='1.2',
        help='max fringe rate for filter')
o.add_option('--inttime',type='float', default='42.9499',
        help='set the inttime specific for lst binned data')
o.add_option('--output',action='store_true',
        help='output each frf as a pkl that can be plotted by plot_frf')
o.add_option('--boxcar',action='store_true',
        help='compute boxcar FIR to compare with frf')
o.add_option('--teff',type='float',
        help='supply tefff for boxcar manually')
o.add_option('--fringe_rate_centered',action='store_true',dest='frc',
        help='center the boxcar around Fringe Rate of maximum optimal FRP')
o.add_option('--noise', '-n', type=float,
        help='noise in Jy to simulate filtering')
o.add_option('--baseline',type=float,
        help='provide integer baseline number to create FRP')
o.add_option('--data_inttime',type=float,default=42.9499,
        help='inttime of data used to recale timebins if inttime option given')
opts,args = o.parse_args(sys.argv[1:])


def get_colors(N):
    '''Returns function with N unique colors'''
    norm=colors.Normalize(vmin=0,vmax=N-1)
    scal_map=cmx.ScalarMappable(norm=norm,cmap='hsv')
    def map_index_to_rgb(index):
        return scal_map.to_rgba(index)
    return map_index_to_rgb
def noise_equivalent_bandwidth(t,H):
    #print 'NEQB amp:', 1./n.sum(n.abs(H)*n.diff(t)[0])**2 
    return 1./n.max(n.abs(H))**2 * n.sum(n.abs(H)**2)*n.diff(t)[0]
    #return 1./n.sum(n.abs(H)*n.diff(t)[0])**2 * n.sum(n.abs(H)**2*n.diff(t)[0])
    #return 1./n.sum(n.abs(H)*n.diff(t)[0])**2 * n.sum(n.abs(H)**2*n.diff(t)[0])

freqs = n.linspace(0.1,0.2,num=203)
aa = a.cal.get_aa(opts.cal, freqs)
nchan = len(freqs)
#pol = a.miriad.pol2str[uv['pol']]

#Get only the antennas of interest
sep2ij, blconj, bl2sep = zsa.grid2ij(aa.ant_layout)

#print "Looking for baselines matching ", opts.ant
#ants = [ b[0] for b in a.scripting.parse_ants(opts.ant, nants) ]
#seps = [ bl2sep[b] for b in ants ]
#seps = n.unique(seps)
PLOT=opts.plot
seps = opts.seps.split(';')
bl_scales = [ float(x) for x in opts.bl_scale.split(',')]
fr_width_scales = [ float(x) for x in opts.fr_width_scale.split(',')]
n_bls=len(bl_scales)
n_width=len(fr_width_scales)
cmap= get_colors(int(1.5*n_bls*n_width))
if opts.chan: mychan = opts.chan
else: mychan = n.floor(nchan/2)

##use calculated inttime to generate correct frf bins
frbins = fringe.gen_frbins(opts.inttime)
#frbins =n.fft.fftshift( n.fft.fftfreq(609,d=opts.inttime))

print 'These are the separations that we are going to use ', seps
print "calculating fringe profile at channel ",mychan
#Get the fir filters for the separation used.
fig_firs,ax_firs=p.subplots(1)
fig_frp,ax_frp=p.subplots(1)
frp={}
for cnt,scale in enumerate(bl_scales):
    firs = {}
    frps = {}
    frp[scale] = {}
    for cn1,frw_scale in enumerate(fr_width_scales):
        frp[scale][frw_scale] ={}
        for sep in seps:
            if opts.baseline:  ij =  a.miriad.bl2ij(opts.baseline); bl =opts.baseline
            else:
                ij_array =  sep2ij[sep].split(',')
                while True:
                    ij = map( int, ij_array.pop().split('_') )
                    bl = a.miriad.ij2bl(*ij )
                    if not blconj[bl]: break
            print "bl_scale = ", scale,"fr_width:", frw_scale, "sep: ",sep,'bl:',bl, 'ant:', ij
            frp[scale][frw_scale][sep], bins = fringe.aa_to_fr_profile(aa, ij, mychan,bins=frbins, pol=opts.pol, bl_scale=scale)
            timebins, firs[sep] = fringe.frp_to_firs(frp[scale][frw_scale][sep], bins, aa.get_afreqs(), fq0=aa.get_afreqs()[mychan], limit_xtalk=True, bl_scale = scale, fr_width_scale = frw_scale, maxfr=opts.maxfr)
            timebins*= opts.data_inttime/opts.inttime

            if False and scale ==1:
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

            frps[sep], frp_freqs = fringe.fir_to_frp(firs[sep],tbins=timebins)
            baselines = ''.join(sep2ij[sep] for sep in seps)

            if PLOT:
                ax_frp.plot(frp_freqs*1e3,frps[sep][mychan]/n.max(frps[sep][mychan]), 
                                label='{0}'.format(scale),color=cmap(cnt))
                ax_firs.plot(timebins,n.abs(firs[sep][mychan]),'-',
                                label='{0}'.format(scale),color=cmap(cnt),
                                linewidth=2)
                ax_firs.plot(timebins,firs[sep][mychan].real,'--',
                                color=cmap(cnt),alpha=.5)
                ax_firs.plot(timebins,firs[sep][mychan].imag,'-.',
                                color=cmap(cnt),alpha=.5)
                ax_firs.set_xlabel('s')
            envelope = n.abs(firs[sep][mychan])
            envelope /= n.max(envelope)
            dt = n.sqrt(n.sum(envelope*timebins**2)/n.sum(envelope))
            dt_50 = (timebins[envelope>0.5].max() - timebins[envelope>0.5].min())
            NEQB_frp = 1.2/noise_equivalent_bandwidth(frp_freqs,frps[sep][mychan])
            #NEQB_fir = noise_equivalent_bandwidth(timebins,firs[sep][mychan])
            print "\t\t\tvariance width [s]:",int(n.round(dt)),
            print "\n\t\t\t50% width = ",int(n.round(dt_50)),
            print "\n\t\t\tT_NEQW (frp) = ",n.round(NEQB_frp)

            if opts.boxcar:
                top_hat = n.zeros_like(firs[sep])
                l_hat = len(top_hat[mychan])
                if opts.teff: box_time = opts.teff
                else: box_time = NEQB_frp
                start = n.round(l_hat/2. - box_time/opts.inttime/2.)
                end = n.round(l_hat/2. + box_time/opts.inttime/2.)
                diff = n.round(box_time/opts.inttime -( end - start))
                if diff != 0: end += diff
                if  (end-start) % 2 == 0: end +=1
                top_hat[:,start:end] += 1.
                t_frp, _ = fringe.fir_to_frp(top_hat,tbins=timebins)
                if opts.frc:
                    t_frp = n.array([ n.roll(t_frp[ch], frps[sep][ch].argmax() - t_frp[ch].argmax()  ,axis=0) for ch in xrange(n.shape(t_frp)[0])] )
                #    t_frp[:,n.logical_or(frp_freqs <= .0002, frp_freqs*1e3 >=1.5)] = 0
                    top_hat = fringe.frp_to_fir(t_frp)
                top_hat /= n.sqrt(n.sum(abs(top_hat)**2,axis=-1).reshape(-1,1))

                NEQB_sinc= 1./noise_equivalent_bandwidth(frp_freqs,t_frp[mychan])

                NEQB_scale = NEQB_sinc/NEQB_frp

                print "\t\t\tBoxcar Len=", (end-start)*opts.inttime,
                print "\n\t\t\tNEQB Boxcar=", n.round(NEQB_sinc,decimals=0)
                print "\t\t\tNEQB scaling=", n.round(NEQB_scale,decimals=4),n.round(1./NEQB_scale,decimals=4)


                if opts.noise:
                    print 'SIMULATING WHITE NOISE at {0}Jy'.format(opts.noise)
                    no_filt_noise = noise((203,711))  * opts.noise
                    no_filt_noise_back = n.copy(no_filt_noise)
                    #wij = n.transpose(f[days[0]][bls_master[0]], [1,0]) #flags (time,freq)
                    wij = gen_flags((203,711),n.zeros((203,711))) #flags (time,freq)
                    #eor1 = noise((21,1451)) * INJECT_SIG #shape of full data
                    dij,wij =n.copy(no_filt_noise) ,n.logical_not(wij)

                    _d = n.convolve( dij[mychan], n.conj(firs[sep][mychan]),mode='same')
                    _w = n.convolve( wij[mychan], n.abs(firs[sep][mychan]),mode='same')
                    _d_frf = n.where(_w > 0, _d/_w,0).reshape((1,len(_d)))

                    _d_bx = n.convolve( dij[mychan], n.conj(top_hat[mychan]),mode='same')
                    _w_bx = n.convolve( wij[mychan], n.abs(top_hat[mychan]),mode='same')
                    _d_bx = n.where(_w_bx>0, _d_bx/_w_bx,0).reshape((1,len(_d)))

                    noise_frf = clip_array(_d_frf,711,axis=1)[0]
                    noise_box = clip_array(_d_bx,711,axis=1)[0]
                    no_filt_noise = clip_array(no_filt_noise,711,axis=1)[mychan]

                    p.figure()
                    p.subplot(311); p.plot(no_filt_noise.real,'k-', alpha=.5)
                    p.subplot(311); p.plot(noise_frf.real, 'k-', label='FRF')
                    p.subplot(311); p.plot(noise_box.real, 'r-', label='Box')
                    p.ylim([-2,2])
                    p.ylabel('Real')
                    p.subplot(312); p.plot(no_filt_noise.imag,'k-', alpha=.5)
                    p.subplot(312); p.plot(noise_frf.imag, 'k-', label='FRF')
                    p.subplot(312); p.plot(noise_box.imag, 'r-', label='Box')
                    p.ylim([-2,2])
                    p.ylabel('Imag')
                    p.subplot(313); p.plot(n.angle(no_filt_noise),'k-', alpha=.5)
                    p.subplot(313); p.plot(n.angle(noise_frf), 'k-', label='FRF')
                    p.subplot(313); p.plot(n.angle(noise_box), 'r-', label='Box')
                    p.ylabel('Phase')

                    p.figure()
                    #p.plot(n.fft.fftshift(n.fft.fft(no_filt_noise.real)),'k-',alpha=.5)
                    p.plot(abs(n.fft.fftshift(n.fft.ifft(noise_box))),'r-',alpha=.8)
                    p.plot(abs(n.fft.fftshift(n.fft.ifft(noise_frf))),'k-',alpha=.8)
                    p.yscale('log')


                if PLOT:
                    if opts.teff: ax_firs.plot(timebins,abs(top_hat[mychan]),'k--',label='T_eff= {0:.0f}'.format(opts.teff))
                    else: ax_firs.plot(timebins,abs(top_hat[mychan]),'k--',label='T_eff= {0:.0f}'.format(NEQB_frp))
                    ax_frp.plot(frp_freqs*1e3,t_frp[mychan]/n.max(t_frp[mychan]),'k--')
                    ax_frp.plot(frp_freqs*1e3,abs(t_frp[mychan]/n.max(t_frp[mychan])),'k--',alpha=.5)
                    #ax_frp.plot(frp_freqs*1e3,t_frp[40]/n.max(t_frp[40]),'g--')
        if opts.output: 
            filename = 'fringe_rate_profile_blscale{scale}_frwidthscale{frwscale}_int{inttime}.pkl'.format(scale=scale,inttime=opts.inttime,frwscale=frw_scale)
            print("saving to: {filename}\n".format(filename=filename))
            f = open(filename,'w')
            out={}
            out['frps'] = frps[sep]
            out['freqs'] = frp_freqs
            out['frp'] = frp[scale][frw_scale][sep]
            out['firs'] = firs[sep]
            out['timebins'] = timebins
            out['chan'] = mychan
            pickle.dump(out,f)
            f.close()
if PLOT:
    for scale in bl_scales:
        for frw_scale in fr_width_scales:
            for sep in seps:
                ax_frp.plot(frp_freqs*1e3,frp[scale][frw_scale][sep],label='frp0 '+sep+' scale:{0}'.format(scale),color='black')
    #ax_frp.plot(frp_freqs*1e3,skew([mn,sq,sk],bins),'k--',label='data')
    ax_frp.set_title('Fitted Fringe Rate Profile')
    ax_frp.set_xlabel('Fringe Rate [mili Hz]')
    #ax_frp.set_xlim([-.7,1.5])
    ax_frp.set_xlim([-.7,5])
    ax_frp.set_ylim([0,1])
    ax_frp.legend(loc='best')

    ax_firs.set_title('Fringe Rates Filter Widths')
    ax_firs.legend(loc='best')
    p.show()

