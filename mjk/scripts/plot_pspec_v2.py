#! /usr/bin/env python
"""
plot the output of pspec_plot_pk_k3pk.py

The plotter also does the final calibration. This just makes a slightly nicer plot using the final calibrated output
npz file.
Pair with plot_p3k_vs_z for a full picture of the power spectrum
"""
def r(number):
    return '%5.0f'%n.round(number,-2)
import matplotlib as mpl
mpl.rcParams['font.size'] = 18
mpl.rcParams['legend.fontsize'] = 14
from pylab import *
import numpy as n,sys,os,re
from capo.cosmo_units import *
from capo import pspec 
import optparse,glob
from scipy.interpolate import interp2d
from capo import pspec,eor_results,cosmo_units
c = 3e8
B =10
k_par_tsys_range = [0.1,1] #range of k_parrallels in which to compute Tsys
print "using default bandwidth of ",B,"MHz for sensitivity line"
#In [3]: F.files
#Out[3]: ['pk', 'kpl', 'k3pk', 'err', 'k3err']
o = optparse.OptionParser()
o.set_usage('plot_pspec.py [options]')
o.set_description(__doc__)
o.add_option('--time',default=0,type=float,
    help='integration time in hours')
opts,args = o.parse_args(sys.argv[1:])

files = sort(args)

#form up the frequency
freqs = []
print "parsing npz file frequencies"
for filename in files:
    print filename,
    try:
        print "npz..",
        freqs.append(n.mean(n.load(filename)['afreqs'])*1e3)
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
print "sorting input files"
files = files[n.argsort(freqs)]
freqs = n.sort(freqs)
print "found freqs"
freqs = n.array(freqs)
print freqs

z = f212z(freqs*1e6)
print "processing redshifts:",z
redshift_files = dict(zip(z,files))
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
if True:
    print "using 50%% reionization Lidz model (%s) and scaling by z"%('power_21cm_z7.32.dat',)
    d = n.array([map(float, L.split()) for L in open('lidz_mcquinn_k3pk/power_21cm_z7.32.dat').readlines()])
    ks, pk = d[:,0], d[:,1]
    k3pk = ks**3 / (2*n.pi**2) * pk     
    print ks.shape,k3pk.shape
    EoR_MODEL_Delta2 = lambda k,z: interp(k,ks,k3pk * mean_temp(z)**2)
    EoR_MODEL_Delta2 = n.vectorize(EoR_MODEL_Delta2)
    if False:
        plot(ks,k3pk*mean_temp(7.83)**2,label='model')
        plot(ks,EoR_MODEL_Delta2(ks,7.8),':',label='my fit')
        legend()
        show();sys.exit()



#LOAD POBER NOISE MODEL
print "loading a few noise model redshifts and building a 2d (z,k) cubic fit model"
# from jcp April 30
#noisefiles = glob.glob('paper32_dcj_0.???.npz')#glob should load them in increasing freq order
#re_f = re.compile(r'paper32_dcj_(\d+\.\d+)\.npz')

# noise from jcp Jul 21
#noisefiles = glob.glob('paper_pessv3_0.???.npz')#glob should load them in increasing freq order
#re_f = re.compile(r'paper_pessv3_(\d+\.\d+)\.npz')

#noisefiles = glob.glob('paper_dcj_lstcnt_pess_0.???.npz')#glob should load them in increasing freq order
#re_f = re.compile(r'paper_dcj_lstcnt_pess_(\d+\.\d+)\.npz')

noisefiles = glob.glob('paper_dcj_lstcnt_sa_pess_0.???.npz')#glob should load them in increasing freq order
re_f = re.compile(r'paper_dcj_lstcnt_sa_pess_(\d+\.\d+)\.npz')

noises = []
noise_ks = []
noise_freqs = []
nk_grid = n.linspace(0,1)*0.5+0.01
for noisefile in noisefiles:
    noise = n.load(noisefile)['T_errs']
    noise_k = n.load(noisefile)['ks']
    bad = n.logical_or(n.isinf(noise),n.isnan(noise))
    noise = noise[n.logical_not(bad)] 
    noise_k = noise_k[n.logical_not(bad)]
    #keep only the points that fall in our desired k range
    noise = noise[noise_k<nk_grid.max()]
    noise_k = noise_k[noise_k<nk_grid.max()]
    print noisefile,n.max(noise),
    noise = n.poly1d(n.polyfit(noise_k,noise,3))(nk_grid)
    noises.append(noise)
    noise_ks.append(nk_grid)
    f = float(re_f.match(noisefile).groups()[0])*1e3 #sensitivity freq in MHz
    print f
    noise_freqs.append(f) 
noise_k_range = [n.min(n.concatenate(noise_ks)),n.max(n.concatenate(noise_ks))]
nk_count = n.mean([len(myks) for myks in noise_ks])
nks = n.linspace(noise_k_range[0],noise_k_range[1],num=nk_count)
noise_interp = n.array([interp(nks,noise_ks[i],noises[i]) for i in range(len(noises))])
NK,NF = n.meshgrid(nks,noise_freqs)
#noise_freqs = n.array(noise_freqs)
POBER_NOISE = interp2d(NK,NF,noise_interp,kind='cubic')#build 2d interpolation model
#FINISHED LOADING POBER NOISE MODEL
if False:
    print noisefiles[3]
    semilogy(noise_ks[3],noises[3],'.',label=noisefiles[3])
    semilogy(noise_ks[3],POBER_NOISE(noise_ks[3],noise_freqs[3]),label='fit')
    grid()
    ylim(1,1e9)
    xlim(0,0.5)
    legend()    
    show()
    sys.exit()

if False:
    ## PRINT A TABLE OF ALL VALUES
    print "%%%"*10
    print "%    TABLE     "
    print "% columns = k [hMpc^-1, k3Pk, k3err, Pober k3Err"
    for i,redshift in enumerate(z):
        print " z = {redshift}".format(redshift=redshift)
        for j in xrange(len(k3Pk[i])):
            print kmags[i][j],
            print k3Pk[i][j],
            print k3err[i][j],
            print 2*POBER_NOISE(kmags [i][j],freqs[i]).squeeze()
    
    print "%%%"*10

#make a nice big 2x4 plot of all the new pspectra
# with Pk on top
# and Delta^2 on the bottom
figure(10,figsize=(17,10))
clf()
Nzs = len(z)-1#all points except P14
print Nzs
j = 0
for i,redshift in enumerate(z):#gotta index z to get the right data
    if n.abs(redshift-7.67)<0.01:continue #skip the P14 point
    k_horizon = n.sqrt(cosmo_units.eta2kparr(30./3e8,redshift)**2 + \
                    cosmo_units.u2kperp(15*freqs[i]*1e6/3e8,redshift)**2)
    print "cutting %d points below k= %4.3f"%(n.sum(kmags[i]<k_horizon),k_horizon)
    window_points = n.argwhere(kmags[i]>k_horizon).squeeze()

    j += 1
    ax_pk = subplot(2,Nzs,j)
    ax_delta = subplot(2,Nzs,Nzs+j)
    #plot PAPER
    #    PAPER Pk

    ax_pk.errorbar(kpars[i],n.abs(Pks[i]),yerr=Pkerr[i],fmt='k+')
    #   PAPER Delta^2
    ax_delta.errorbar(kmags[i][window_points],k3Pk[i][window_points],yerr=k3err[i][window_points],fmt='k.', capsize=0)

    #plot PAPER sensitivity
    ax_delta.plot(kmags[i],POBER_NOISE(kmags[i],freqs[i]),'--k')


    #plot MWA 32T data
    #print "loading MWA data near redshift",redshift
    MWA_z,MWA_pspec = eor_results.z_slice(redshift,eor_results.MWA_32T_all())
    #print "found z = ",MWA_z
    if n.abs(MWA_z-redshift)<0.5:
        ax_delta.plot(MWA_pspec[:,0],MWA_pspec[:,2],'kd')
    #plot the GMRT paciga 2014 data
    GMRT_z,GMRT_pspec = eor_results.z_slice(redshift,eor_results.GMRT_2014_all())
    #print "loading GMRT 2014 data near redshift:",GMRT_z
    if n.abs(GMRT_z - redshift)<0.5:
        ax_delta.plot(GMRT_pspec[:,0],GMRT_pspec[:,2],'kx')

    #plot the lidzian model
    ax_delta.plot(kmags[i],EoR_MODEL_Delta2(kmags[i],redshift),'k')

    #ax_pk.text(-0.4,8e9,"z = %3.2f"%redshift)
    ax_pk.set_title("z = %3.2f"%redshift,size=14)
    #setup the axes and such
    # set up PK axis
    ax_pk.set_yscale('log',nonposy='clip')
    if i==0:ax_pk.set_ylabel('P($k_\perp$=%3.2f) [mk$^2$/Mpc$^3$]'%kperps[i])
    else: ax_pk.set_yticklabels([])
    ax_pk.set_xlabel('$k_\parallel$ [h Mpc$^{-1}$]')
    ax_pk.vlines([-cosmo_units.eta2kparr(30./3e8,redshift),cosmo_units.eta2kparr(30/3e8,redshift)],1e6,1e17,linestyles='--')
    ax_pk.set_ylim([1e5,1e10])
    ax_pk.set_xlim([-0.6,0.6])
    ax_pk.grid()
    # set up Delta axis
    ax_delta.set_yscale('log',nonposy='clip')
    if i==0:ax_delta.set_ylabel('$k^3 / 2 \pi^2 P(k)$ [mK$^2$]')
    else: ax_delta.set_yticklabels([])
    ax_delta.set_xlabel('k [h Mpc$^{-1}$]')
    ax_delta.set_ylim([1,1e8])
    ax_delta.set_xlim([0,0.6])
    ax_delta.grid()
tight_layout()
subplots_adjust(wspace=0)
draw()
savefig('PAPER_32T_new_redshifts.png')
#plot some k bins versus redshift
ratio = []
sigmas = []
for i,k in enumerate([0.2,0.3,0.4]):
    noise_vs_z = POBER_NOISE(k,freqs).squeeze()
    myk_slice_indexes =     [n.abs(k-kmags[j]).argmin() for j in range(len(kmags))]
    k3Pk_err_k_slice =      [k3err[j][myk_slice_indexes[j]] for j in range(len(kmags))]
    k3Pk_k_slice =          [k3Pk[j][myk_slice_indexes[j]] for j in range(len(kmags))]
    for j,redshift in enumerate(z):
        #ratio.append(n.sqrt(2)*noise_vs_[j]/k3Pk_err_k_slice[j])
        #print k,redshift,"noise theoreti0cal/bootstrap = ",ratio[-1]
        sigmas.append(k3Pk_k_slice[j]/k3Pk_err_k_slice[j])
        print k,"nsigmas = ",sigmas[-1]
#print n.mean(ratio)

print n.sqrt(n.sum(n.array(sigmas)**2))
#TODO What are the units of the noise!?!?! The error bars I get 

figure(4)
for i,k in enumerate([0.1,0.2,0.3,0.4]):
    clf()
    axz=subplot(111)
    axz.set_yscale('log',nonposy='clip')
    axf = axz.twiny() #switch to the upper 
    #plot the sensitivity
    noise_vs_z = POBER_NOISE(k,freqs).squeeze()
    noise_line = plot(freqs,noise_vs_z,'--k')[0]
    xlabel('frequency [MHz]')


    #plot the model
    #plot LIDZIAN model
    plot(freqs,EoR_MODEL_Delta2(k,z[::-1])[::-1],'k')




    #find the k points closest to the desired k vs redshift
    myk_slice_indexes =     [n.abs(k-kmags[j]).argmin() for j in range(len(kmags))]
    k3Pk_k_slice =          [k3Pk[j][myk_slice_indexes[j]] for j in range(len(kmags))]
    k3Pk_err_k_slice =      [k3err[j][myk_slice_indexes[j]] for j in range(len(kmags))]
    myk_slice =             [kmags[j][myk_slice_indexes[j]] for j in range(len(kmags))]
    #print "k slice",k,"hMpc^-1"
    for j,redshift in enumerate(z):
        if n.abs(redshift-7.67)<0.01: #plot the arp point special
            errorbar(freqs[j],k3Pk_k_slice[j],yerr=k3Pk_err_k_slice[j],fmt='.k',capsize=0)
        else:
            errorbar(freqs[j],k3Pk_k_slice[j],yerr=k3Pk_err_k_slice[j],fmt='db',capsize=0,mec='b')#plot the PAPER data points
            print k,'&',
            print n.round(redshift,2),'&',
            print r(k3Pk_k_slice[j]),'&',
            print '$\pm$',r(n.sqrt(2)*noise_vs_z[j]),'&',
            print n.round(k3Pk_k_slice[j]/(n.sqrt(2)*noise_vs_z[j]),1),'&',
            if (n.abs(k3Pk_k_slice[j])-n.sqrt(2)*noise_vs_z[j])<0:
                print "ULim",#print "theory error consistent with 0 at 95% conf"
            else:
                print "Det",#print "theory error not consistent with 0 at 95% conf"
            print '&',
            print '$\pm$',r(k3Pk_err_k_slice[j]),'&',
            print n.round(k3Pk_k_slice[j]/k3Pk_err_k_slice[j],1),'&',
            if (n.abs(k3Pk_k_slice[j]) - k3Pk_err_k_slice[j])<0: 
                print "ULim", # print "boostrap error consistent with 0 at 95% conf",
            else:
                print "Det",#print "bootstrap error not constistent with 0 at 95% conf"
            print '\\tabularnewline'
        
    #plot the MWA data
    MWAz,MWAp = eor_results.MWA_32T_at_k(k)
    #plot the GMRT paciga 2014 data
    GMRT_z,GMRT_pspec = eor_results.k_slice(k,eor_results.GMRT_2014_all())
    plot(1421/(GMRT_z+1),GMRT_pspec[:,1],'+k')
    ylim([1,5e4])

    ylabel('$k^3/2\pi^2 P(k=%3.2f)$  [mK$^2$]'%(k))
#    sca(axz)
    #put matching redshifts on the bottom of the axis
    
    freq_ticks = axf.get_xticks()#n.round(1421./(1+z_ticks)).astype(int)
    z_ticks = n.round(1421/freq_ticks -1) 
    axz.set_xlim(axf.get_xlim())#the underlying axis is redshift, the labels are freqs
    axz.set_xticks(1421/(z_ticks+1))
    axz.set_xticklabels(z_ticks)
    axz.set_xlabel('redshift')


    axz.grid(which='both')
    #print 'psa32_pspec_k_%4.2f_log.png'%k
    savefig('psa32_pspec_k_%4.2f_log.png'%k)
    #change to linear units
    noise_line.set_visible(False)#turn off dashed line
    axf.set_yscale('linear',nonposy='clip')
    axf.set_ylim([-2*n.max(k3Pk_k_slice),2*n.max(k3Pk_k_slice)])
    axz.set_yscale('linear',nonposy='clip')
    axz.set_ylim([-2*n.max(k3Pk_k_slice),2*n.max(k3Pk_k_slice)])
#    axz.set_ylim([-2*n.max(k3Pk_k_slice),2*n.max(k3Pk_k_slice)])

    #ylim([-5e4,5e4])
    fill_between(freqs,noise_vs_z,y2=-1*noise_vs_z,color='k',alpha=0.3,zorder=-10)
    #print 'psa32_pspec_k_%4.2f_lin.png'%k
    axz.set_ylabel('$k^3/2\pi^2 P(k=%3.2f)$  [mK$^2$]'%(k))
    draw()
    savefig('psa32_pspec_k_%4.2f_lin.png'%k)
#plot the power spectrum at each redshift
for i,redshift in enumerate(z):
    #THE LINEAR FIGURE
    figure(1)
    clf()
    ax2 = subplot(111)
    #plot new PAPER data
    errorbar(kmags[i],k3Pk[i],yerr=k3err[i],fmt='k.', capsize=0)
    fill_between(kmags[i],POBER_NOISE(kmags[i],freqs[i]),y2=-POBER_NOISE(kmags[i],freqs[i]),color='k',alpha=0.3)
    ylim([-1e5,1e5])
    xlim([0,0.6])
    #plot the predicted noise level
    Tsys = (450*(freqs[i]/160.)**(-2.6) + 100)*1e3
    print "Tsys = ",Tsys
    print "B = ",B

    xlabel('k [h Mpc$^{-1}$]')
    ylabel('$k^3 / 2 \pi^2 P(k)$ [mK$^2$]')
    grid()
    tight_layout()
    print 'pspec_lin_z_%4.2f.png'%redshift
    savefig('pspec_lin_z_%4.2f.png'%redshift)


    #the LOG FIGURE
    figure(2)
    clf()
    ax1 = subplot(121) #plot Pk
    ax1.set_yscale('log',nonposy='clip') #note, the nonposy prevents error bars from disappearing if negative
    #plot +s over the fg points from arps previous paper
    fg_points = n.argwhere(Pks[i]>1e14).squeeze()
    if len(fg_points)>0:
        errorbar(kpars[i][fg_points],n.abs(Pks[i][fg_points]),yerr=Pkerr[i][fg_points],fmt='kx')     
    errorbar(kpars[i],n.abs(Pks[i]),yerr=Pkerr[i],fmt='k.') #plot the PAPER data
    #plot(kpars[i],n.ones_like(kpars[i])*noise_Pk,':k') #plot the noise
    ylim([1e6,1e17])  
    ylabel('P($k_\perp$=%3.2f) [mk$^2$/Mpc$^3$]'%kperps[i])
    xlabel('$k_\parallel$ [h Mpc$^{-1}$]')
    print "sky delay @ redshift",redshift,"=",cosmo_units.eta2kparr(30./3e8,redshift)
    vlines([-cosmo_units.eta2kparr(30./3e8,redshift),cosmo_units.eta2kparr(30/3e8,redshift)],1e6,1e17,linestyles='--')
    grid()
    text(-0.4,1e15,"z = %3.2f"%redshift)
    ax2 = subplot(122) #plot Delta^2
    ax2.set_yscale('log',nonposy='clip')
    #plot MWA 32T data
    #print "loading MWA data near redshift",redshift
    MWA_z,MWA_pspec = eor_results.z_slice(redshift,eor_results.MWA_32T_all())
    #print "found z = ",MWA_z
    if n.abs(MWA_z-redshift)<0.5:
        plot(MWA_pspec[:,0],MWA_pspec[:,2],'kd')
    #plot the GMRT paciga 2014 data
    GMRT_z,GMRT_pspec = eor_results.z_slice(redshift,eor_results.GMRT_2014_all())
    #print "loading GMRT 2014 data near redshift:",GMRT_z
    if n.abs(GMRT_z - redshift)<0.5:
        plot(GMRT_pspec[:,0],GMRT_pspec[:,2],'kx')


    #plot PAPER data
    errorbar(kmags[i],k3Pk[i],yerr=k3err[i],fmt='k.', capsize=0)
    #plot PAPER noise
    #plot(kmags[i],noise_pspec,':k')
    #plot(POBER_NOISE_ks,POBER_NOISE_Delta2,'--')#plot the noise from Jonnie
    plot(kmags[i],POBER_NOISE(kmags[i],freqs[i]),'--k')

    #plot LIDZIAN model
    plot(kmags[i],EoR_MODEL_Delta2(kmags[i],redshift),'k')
    #print "peak of the Lidzian model", n.max(EoR_MODEL_Delta2(kmags[i],redshift))
    xlabel('k [h Mpc$^{-1}$]')
    ylabel('$k^3 / 2 \pi^2 P(k)$ [mK$^2$]')
    ylim([1,1e8])
    xlim([0,0.6])
    grid()
    tight_layout()
    if False:
        figure(3)
        semilogy(kpars[i],Pkerr[i],'.k')
        ylabel('2sigma error on P($k_\perp$=%3.2f) [mk$^2$/Mpc$^3$]'%kperps[i])
        print "kperp = ", kperps[i]
        xlabel('$k_\parallel$ [h Mpc$^{-1}$]')
        ylim([1e6,5e7])
        grid(which='both')
    print 'pspec_log_z_%4.2f.png'%redshift
    savefig('pspec_log_z_%4.2f.png'%redshift)
    
#show()

