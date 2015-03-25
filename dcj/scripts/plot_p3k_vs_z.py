#! /usr/bin/env python
"""
Plot a z vs kpl waterfall using the npz files output by pspec_plot_k3pk.py
plot_pk_waterfall.py *npz
"""
import matplotlib as mpl
mpl.rcParams['font.size'] = 18
mpl.rcParams['legend.fontsize'] = 14
from pylab import *
import numpy as n,sys,os,re
from capo.cosmo_units import *
from capo import pspec 
import optparse
c = 3e8

#In [3]: F.files
#Out[3]: ['pk', 'kpl', 'k3pk', 'err', 'k3err']
o = optparse.OptionParser()
o.set_usage('plot_p3k_vs_z.py [options]')
o.set_description(__doc__)
o.add_option('-k',default=0.2,type=float,
    help='The input k at which to plot. [hMpc^-1, default=0.2]')
opts,args = o.parse_args(sys.argv[1:])
myk = opts.k

files = sort(args)

#form up the frequency
dchan = n.array([int(F.split('/')[1].split('_')[1])-int(F.split('/')[1].split('_')[0]) for F in files])
chans = n.array([int(F.split('/')[1].split('_')[0]) for F in files])+dchan/2
I = n.argsort(chans)
chans = chans[I]
freqs = chans/2. + 100
print freqs
z = f212z(freqs*1e6)
umags = 30/(c/freqs*1e6)
kperps = umags*pspec.dk_du(z)
#k = n.sqrt(kpl**2 + kpr**2)
#k3 = n.abs(k**3 / (2*n.pi**2))
#print [len(p) for p in P]


#plot a single k bin
#mykpl = 0.35
#assume 30m baselines at 150Mhz
#print "plotting k_parallel = ",mykpl
Pkk = []
Pkk_err = []
k3Pk = []
k3err = []
neg = []
for i,FILE in enumerate(files):
    F = n.load(FILE)
#    print FILE,z[i]
    #ki.append(n.abs(F['kpl'] - mykpl).argmin())
    k = n.sqrt(F['kpl']**2 + kperps[i]**2)
    ki = n.abs(k-myk).argmin()
    ki_neg = n.abs(k+myk).argmin()
#    print ki,ki_neg,'/',len(k)
    #P = (F['k3pk'][ki] + F['k3pk'][ki_neg])/2
    #Perr =n.sqrt((F['k3err'][ki]**2 + F['k3err'][ki_neg]**2)/2)
    P = F['k3pk'][ki]
    Perr =F['k3err'][ki]
    neg.append(0)
    if P<0:
        if (P-Perr)<0:
            neg[-1] = 1
        P = n.abs(P)
    k3Pk.append(P)
    k3err.append(Perr)
    #print F['kpl'][ki],P,Perr
    #kineg.append(n.abs(F['kpl'] + mykpl).argmin())
    #Pkk.append((F['pk'][ki] + F['pk'][kineg])/2)
    #Pkk_err.append(n.sqrt((F['pk'][ki]**2 + F['pk'][kineg]**2)/2))
#Pkk = n.array(Pkk)
#Pkk_err = n.array(Pkk_err)
myk = k[ki]
k3Pk = n.array(k3Pk)
k3err = n.array(k3err)
#Pk_at_mykpl =      n.array( [P[i,ki[i]] for i in range(len(ki))])
#Pk_err_kpl =       n.array(  [Perr[i,ki[i]] for i in range(len(ki))])
#Pk_at_mykpl_neg =  n.array(  [P[i,kineg[i]] for i in range(len(kineg))])
#Pk_err_kpl_neg =   n.array(  [Perr[i,kineg[i]] for i in range(len(kineg))])



#Load Lidzian models
def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

import glob
re_z = re.compile(r'power_21cm_z(\d+\.\d+)\.dat')
model_slice =  []
for filename in glob.glob('lidz_mcquinn_k3pk/power*dat'):
#    print 'Reading', filename
    d = n.array([map(float, L.split()) for L in open(filename).readlines()])
    ks, pk = d[:,0], d[:,1]
    z_file = float(re_z.match(os.path.basename(filename)).groups()[0])+0.5
    umag = 16*(1+z_file)/7.6
    kperp =  umag * pspec.dk_du(z_file)  
    #k = n.sqrt(kperp**2 + mykpl**2)
#    z = C.pspec.f2z(.160)
    k3pk = ks**3 / (2*n.pi**2) * pk * mean_temp(z_file)**2
    ki = n.abs(ks-myk).argmin()
    model_slice.append([z_file,k3pk[ki]])
model_slice = n.array(model_slice)
I = n.argsort(model_slice[:,0])
model_slice = model_slice[I]


figure(figsize=(10,5))
ax1 = subplot(111)
#ax1.errorbar(freqs,Pk_at_mykpl,yerr=Pk_err_kpl)
#ax1.errorbar(freqs,Pk_at_mykpl_neg,yerr=Pk_err_kpl_neg)
#Pkk = (Pk_at_mykpl+Pk_at_mykpl_neg)/2
#Pkk_err = n.sqrt((Pk_err_kpl**2 + Pk_err_kpl_neg**2)/2)
ax1.set_yscale('log',nonposy='clip')
#ax1.errorbar(freqs,Pkk,yerr=Pkk_err,fmt='xk')
#ax1.errorbar(freqs,n.sqrt(k3Pk),yerr=n.sqrt(k3err),fmt='xb')
arrowlen=10
label='_nolegend_'
islabeled = False
for F,P,dP in zip(freqs,k3Pk,k3err):
    P = n.sqrt(P)
    dP = n.sqrt(dP)
    #annotate('',xy=[F,(P+dP)/arrowlen],xycoords='data',xytext=[F,P+dP],textcoords='data',
    #    arrowprops={'arrowstyle':'->','color':'b','lw':1})
    if F==164.5:
        fmt='k'
        lw=1
        label='Parsons2014'
    else:
        fmt='db'
        lw=3
        if not islabeled:
            label='this work'
            islabeled=True
        else:
            label='_nolegend_'
    if n.isnan(P):
        print "...skipping"
        continue
    ax1.errorbar(F,P,yerr=dP,fmt=fmt,lw=lw,capsize=10,mew=lw,mec='b',ms=15,
            label=label)


#plot the Lidzian Model
ax1.plot(1421/(1+model_slice[:,0]),n.sqrt(model_slice[:,1]),'k',lw=4)
ax1.set_ylim([0.1,2e3])
ax1.set_ylabel('$\sqrt{k^3/2\pi^2 P(k=%3.2f)}$  [mK]'%(myk))

ax1.set_xlim([140,178])
xlims = ax1.get_xlim()
ax2 = ax1.twiny()

x_z = n.array(list(set(n.round(f212z(n.linspace(ax1.get_xlim()[0],ax1.get_xlim()[1])*1e6))))+[7,6])
x_zfreqs = f21/(1+x_z)/1e6
ax2.set_xticks(x_zfreqs)
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticklabels(x_z)
ax1.set_xlabel('freq [MHz]')
ax2.set_xlabel('redshift')
#A crude estimate of MWA sensitivity at k=0.2 from Beardsley et al 2012
freqs = n.linspace(110,200)
z = f212z(freqs*1e6) 
dT_MWA = 20e-3*28.*((1+z)/10)**(1/2)*(freqs/150.)**-2.6 #mK^2
ax1.plot(freqs,n.sqrt(dT_MWA),'c',label='MWA sensitivity')

#A crude estimate of PAPER sensitivity at k=0.2 from Parsons et al 2012a
freqs = n.linspace(120,180) 
dT_PAPER = 20 * (freqs/150.)**-2.6 #mK^2
ax1.plot(freqs,n.sqrt(dT_PAPER),'m',label='PAPER sensitivity')

#plot the GMRT results
GMRT_paciga_2013 = n.array([[0.1,2e5],[0.13,4e5],[0.16,1e5],[0.19,1.9e5],[0.275,2e5],[0.31,4e5],[0.4,6e5],[0.5,8e4]])
i_GMRTk = n.abs(GMRT_paciga_2013[:,0]-myk).argmin()
GMRT_dk =  n.abs(GMRT_paciga_2013[:,0]-myk).min() 
if GMRT_dk<0.0125: #only plot GMRT points when they are within 0.0125 hMpc^-1
    print "GMRT i_k,k",i_GMRTk,GMRT_paciga_2013[i_GMRTk,0]
    ax1.errorbar(1421/(1+8.6),n.sqrt(GMRT_paciga_2013[i_GMRTk,1]),yerr=n.sqrt(1.9e5/2),fmt='xk',label='Paciga2013')


ax1.legend(numpoints=1,loc='lower right')
savefig('psa_p3k_vs_z_%d.png'%(myk*100))
