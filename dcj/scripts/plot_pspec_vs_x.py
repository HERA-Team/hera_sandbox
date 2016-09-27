#! /usr/bin/env python
"""
plot the output of pspec_plot_pk_k3pk.py vs a variable

Takes as input a file listing npz pspecs and the relevant variable
ex:
usage: plot_pspec_vs_x.py pspec_vs_frwidth.txt

pspec_vs_frwidth.txt
#FILE frwidth
Jun3_optimal_frwidth10.0/30_50/I/pspec_Jun3_optimal_frwidth10.0_30_50_I.npz   10
Jun3_optimal_frwidth1.48/30_50/I/pspec_Jun3_optimal_frwidth1.48_30_50_I.npz   1.48
Jun3_optimal_frwidth1.8/30_50/I/pspec_Jun3_optimal_frwidth_1.8_30_50_I.npz    1.8
Jun3_optimal_frwidth2.0/30_50/I/pspec_Jun3_optimal_frwidth2.0_30_50_I.npz     2.0
Jun3_optimal_frwidth4.0/30_50/I/pspec_Jun3_optimal_frwidth4.0_30_50_I.npz     4.0

"""
def r(number):
    return '%5.0f'%n.round(number,1)
    #return number

import matplotlib as mpl
mpl.rcParams['font.size'] = 12
mpl.rcParams['xtick.labelsize']= 18
mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['lines.markersize']= 12

mpl.rcParams['legend.numpoints']  = 1
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['figure.dpi'] = 500
mpl.rcParams['savefig.dpi'] = 500
mpl.rcParams['savefig.format'] ='png'
from pylab import *
import numpy as n,sys,os,re
from capo.cosmo_units import *
import optparse,glob
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
from capo import pspec,eor_results,cosmo_units,zsa
import ipdb
import aipy as a
import capo, capo.frf_conv as fringe
import matplotlib.pyplot as p
from IPython import embed

c = 3e8
B =10
k_par_tsys_range = [0.1,1] #range of k_parrallels in which to compute Tsys

def noise_equivalent_bandwidth(t,H):
    return 1./n.max(n.abs(H))**2 * n.sum(n.abs(H)**2)*n.diff(t)[0]

print "using default bandwidth of ",B,"MHz for sensitivity line"
#In [3]: F.files
#Out[3]: ['pk', 'kpl', 'k3pk', 'err', 'k3err']
o = optparse.OptionParser()

o.set_usage('plot_noise_vs_x.py [options]')
a.scripting.add_standard_options(o, cal=True)
o.set_description(__doc__)
o.add_option('--scale_noise',type=float,
    help='if set, scale the sensitivity line by <scale_noise>/x. Off by default')
#o.add_option('--time',default=0,type=float,
#    help='integration time in hours')
#o.add_option('--title',type='string',
#        help='supply title for output')
#o.add_option('--inttime', default=42.9499,type=float,
#    help='inttime to use in making frf filter (default = 42.9499s)')
opts,args = o.parse_args(sys.argv[1:])



#read in the input file
lines = open(args[0]).readlines()
#clean off any empty lines
lines = [line for line in lines if len(line.split())>1]
files = [line.split()[0] for line in lines if not line.startswith('#')]
x = np.array(map(float,[line.split()[1] for line in lines if not line.startswith('#')]))
try:
    print "looking for noise scale factor in third column..",
    noise_scale = np.array(map(float,[line.split()[2] for line in lines if not line.startswith('#')]))
    print "[FOUND]"
except(IndexError):
    print "[NOT FOUND]"
    noise_scale = np.ones_like(x)
print "looking for the variable name.."
if lines[0].startswith('#FILE'):
    varname=lines[0].split()[1]
    print "  found",varname
else:
    varname='x'
    print "  not found, using",varname
    print "  set first line of input file to read #FILE <varname>, to name your variable"
filetitle=args[0][:-4] #strip off the .txt
print "using title = ",filetitle

#get the freqs from the filenames
freqs = []
inttimes = []
chans = []
print "parsing npz file frequencies"
for filename in files:
    chans = filename.split('/')[1]
    print filename,
    try:
        print "npz..",
        freqs.append(n.load(filename)['freq']*1e3)
        print "[Success]"
    except(KeyError):
        print "[FAIL]"
        try:
            print "looking for path like RUNNAME/chan_chan/I/pspec.npz"
            dchan = int(filename.split('/')[1].split('_')[1])-int(filename.split('/')[1].split('_')[0])
            chan = int(filename.split('/')[1].split('_')[0]) + dchan/2
            freqs.append(chan/2. + 100) #a pretty good apprximation of chan 2 freq for 500kHz channels
        except(IndexError):
            print "[FAIL] no freq found. Skipping..."
            continue


print "found ",len(x),"input files"
x = n.array(x)
#print n.sort(x), max_frf_freqs*1e3
print "sorting input files"
print files,np.argsort(x)
files = [files[xi] for xi in n.argsort(x)]
freqs = n.array(freqs)
freqs = freqs[n.argsort(x)]
noise_scale = noise_scale[np.argsort(x)]
print "Power Spectrum at {0:.2f}Mhz".format(freqs[0])

x = n.sort(x)

print "loading files",files
z = f212z(freqs*1e6)
print "processing redshifts:",z
#redshift_files_pads = dict(zip(z,files,x))
umags = 30/(c/(freqs*1e6))
print "umags = ",umags
kperps = umags*pspec.dk_du(z)
Pkk = []
Pkk_err = []
k3Pk = []
k3err = []
neg = []
kpars = []
kmags = []
Pks = []
Pkerr =[]
for i,FILE in enumerate(files):
    F = n.load(FILE)
    print FILE,z[i]
    k = n.sqrt(F['kpl']**2 + kperps[i]**2)
    k3Pk.append(F['k3pk'])
    Pks.append(F['pk'])
    k3err.append(F['k3err'])
    kpars.append(F['kpl'])
    kmags.append(F['k'])
    Pkerr.append(F['err'])
#clean off this stupid extra point
for i in xrange(len(kmags)):
    if n.round(z[i],2)==7.68:continue
    k0 = n.argwhere(kpars[i]==0).squeeze()
    print "deleting dumb point",k0," out of ",len(kpars[i])
    Pks[i] = n.delete(Pks[i]  ,k0,0) 
    kpars[i] = n.delete(kpars[i],k0,0)
    Pkerr[i] = n.delete(Pkerr[i],k0,0)

#Load Lidzian models
def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK



#LOAD POBER NOISE MODEL
#print "loading a few noise model redshifts and building a 2d (z,k) cubic fit model"
# from jcp April 30
#noisefiles = glob.glob('paper32_dcj_0.???.npz')#glob should load them in increasing freq order
#re_f = re.compile(r'paper32_dcj_(\d+\.\d+)\.npz')

# noise from jcp Jul 21
#noisefiles = glob.glob('paper_pessv3_0.???.npz')#glob should load them in increasing freq order
#re_f = re.compile(r'paper_pessv3_(\d+\.\d+)\.npz')

#noisefiles = glob.glob('paper_dcj_lstcnt_pess_0.???.npz')#glob should load them in increasing freq order
#re_f = re.compile(r'paper_dcj_lstcnt_pess_(\d+\.\d+)\.npz')

#noisefiles = glob.glob('paper_dcj_lstcnt_sa_pess_0.???.npz')#glob should load them in increasing freq order
noisefiles = glob.glob('../21cmsense_noise/psa6240_v003drift_mod_0.???.npz')#glob should load them in increasing freq order
noisefiles.sort()
#re_f = re.compile(r'paper_dcj_lstcnt_sa_pess_(\d+\.\d+)\.npz')
re_f = re.compile('../21cmsense_noise/psa6240_v003drift_mod_(\d+\.\d+)\.npz')#glob should load them in increasing freq order


#build a redshift independent noise model
noises = []
noise_ks = []
noise_freqs = []
nk_grid = n.linspace(0,1)*0.5+0.01
for noisefile in noisefiles:
    noise = n.load(noisefile)['T_errs']
    noise_k = n.load(noisefile)['ks']

    #look for infs at low k and re-assign large values
    #noise[ n.logical_and(n.isinf(noise), noise_k < .15)] = 1e+6

    bad = n.logical_or(n.isinf(noise),n.isnan(noise))
    noise = noise[n.logical_not(bad)] 
    noise_k = noise_k[n.logical_not(bad)]
    #keep only the points that fall in our desired k range
    good_k = noise_k < nk_grid.max()
    noise = noise[noise_k<nk_grid.max()]
    noise_k = noise_k[noise_k<nk_grid.max()]
    print noisefile,n.max(noise),

    noise_k = n.insert(noise_k,0,0)
    noise = n.insert(noise,0, n.min([1e+3,n.median(noise)]))
    #noise = n.insert(noise,0, 0)

    #nk3 = noise_k**3/(2*np.pi**2)
    #tmp = n.linalg.lstsq( n.vstack([nk3,n.zeros((3,len(nk3)))]).T,noise)
    tmp = n.polyfit(noise_k,noise,3,full=True)
    noise = n.poly1d(tmp[0])(nk_grid)
    noises.append(noise)
    noise_ks.append(nk_grid)
    f = float(re_f.match(noisefile).groups()[0])*1e3 #sensitivity freq in MHz
    print f
    noise_freqs.append(f) 
    #embed()

noise_k_range = [n.min(n.concatenate(noise_ks)),n.max(n.concatenate(noise_ks))]
nk_count = n.mean([len(myks) for myks in noise_ks])
nks = n.linspace(noise_k_range[0],noise_k_range[1],num=nk_count)
noise_interp = n.array([interp(nks,noise_ks[i],noises[i]) for i in range(len(noises))])
NK,NF = n.meshgrid(nks,noise_freqs)
#noise_freqs = n.array(noise_freqs)
POBER_NOISE = interp2d(NK,NF,noise_interp,kind='linear')#build 2d interpolation model
#FINISHED LOADING POBER NOISE MODEL


#plot a few k bins vs the input parameter,
# save each as a seperate file

Nzs = len(z)
Nint = len(x)
print 'Number of x: {0}'.format(Nint)

for nk,k in enumerate([.2,.3,.4]):
    print '\tPlotting for k={0}'.format(k)
    figure(10,figsize=(17,10))
    clf()
    ax_delta = subplot(1,1,1)
    #Nzs = len(z)-1#all points except P14

    myk_slice_indexes =     [n.abs(k-kmags[j]).argmin() for j in range(len(kmags))]
    k3Pk_err_k_slice =      [k3err[j][myk_slice_indexes[j]] for j in range(len(kmags))]
    k3Pk_k_slice =          [k3Pk[j][myk_slice_indexes[j]] for j in range(len(kmags))]

    #tmp_noise = 2*POBER_NOISE( k , freqs[0])*
    #tmp_noise = POBER_NOISE( k , freqs[0] )
    tmp_noise = 2*POBER_NOISE( k , freqs[0] ) * noise_scale

    #ax_delta.axhline(tmp_noise,color='k',linestyle='--')
    ax_delta.plot(x,tmp_noise,color='k',linestyle='--')


    for i,itime in enumerate(x):#gotta index z to get the right data
        redshift = z[i]
        x_var = x[i]
        #   PAPER Delta^2
        if k3Pk_k_slice[i] > 0:
            ax_delta.errorbar(x_var,k3Pk_k_slice[i],yerr=k3Pk_err_k_slice[i],fmt='.k',capsize=0,mec='k')#plot the PAPER data points
            #  plot the negative points
        else:
            ax_delta.errorbar(x_var,
                            n.abs(k3Pk_k_slice[i]),
                            yerr=k3Pk_err_k_slice[i],fmt='.',color='0.5', capsize=0)


    #ax_pk.text(-0.4,8e9,"z = %3.2f"%redshift)
    ax_delta.set_ylabel('$k^3 / 2 \pi^2 P(k)$ [mK$^2$]')
    #setup the axes and such

    ax_delta.set_title("k = %3.2f [hMpc^-1]"%k,size=14)
    ax_delta.set_ylabel('$k^3 / 2 \pi^2 P(k)$ [mK$^2$]')
    ax_delta.set_yscale('log',nonposy='clip')
    ax_delta.set_xlabel(varname)
    ax_delta.set_xlim([x.min()*.9,x.max()*1.1])
    ax_delta.set_ylim([1,1e5])
    ax_delta.grid()
    tight_layout()
    subplots_adjust(wspace=0)
    draw()

    savefig(filetitle+'_k{0}_'.format(k)+chans+'.png')
