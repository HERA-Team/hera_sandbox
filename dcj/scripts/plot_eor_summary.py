#! /usr/bin/env python
'''
Creates a summary plot of the Best data vs Redshift from  capo.eor_results
'''
import matplotlib as mpl
fontsize=20
textsize=10
#mpl.rcParams['font.size']  = 20
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['legend.numpoints']  = 1
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.fontsize'] = fontsize
mpl.rcParams['figure.dpi'] = 400
mpl.rcParams['savefig.dpi'] = 400
mpl.rcParams['savefig.format'] ='png'
mpl.rcParams['lines.markeredgewidth'] = 0
mpl.rcParams['lines.markersize'] = 9
#mpl.rcParams['lines.markersize'] = 7

capsize = 3.5
plot_psa32 = True
plot_psa32_multiz = True
plot_psa64 = True
plot_psa64_multiz=False
from capo import eor_results,cosmo_units,pspec
import os,sys, numpy as n, matplotlib.pyplot as p
import optparse
import ipdb

c=299792458.

o=optparse.OptionParser()
o.set_usage('plot_eor_summary.py [options]')
o.set_description(__doc__)
o.add_option('--plot', action='store_true',
        help='outputs plot before saving file')
o.add_option('--models', type='str',
        help='a string that globs a list of 21cmfast pspec output files')

opts,args = o.parse_args(sys.argv[1:])

files = args
freqs=[]
for filename in files:
    try:
        freqs.append(n.load(filename)['freq']*1e3)
        print('Frequency found: {0}'.format(freqs[-1]))
    except(KeyError):
        try:
            dchan = int(filename.split('/')[-1].split('.')[0].split('_')[2])-int(filename.split('.')[0].split('_')[1])
            chan = int(filename.split('/')[-1].split('.')[0].split('_')[1]) + dchan/2.
            freqs.append(chan/2. + 100) #a pretty good apprximation of chan 2 freq for 500kHz channels
            print('Frequency found: {0}'.format(freqs[-1]))
        except:
            print('No frequency found for file: '+filename+'\n')
            print('Skipping')
            continue
freqs = n.array(freqs)
zs = pspec.f2z(freqs*1e-3)

fig = p.figure(figsize=(10,5))
ax = fig.add_subplot(111)



#plot the GMRT paciga 2014 data
GMRT = eor_results.GMRT_2014_all()
GMRT_results = {}
print('GMRT')
for i,z in enumerate(GMRT.keys()):
    #index = n.argwhere(GMRT[z][:,0] - .2 < 0.1).squeeze()
    freq= pspec.z2f(z)
    k_horizon = n.sqrt(cosmo_units.eta2kparr(30./c,z)**2 + \
                    cosmo_units.u2kperp(15*freq*1e6/c,z)**2)
    index = n.argwhere( abs(GMRT[z][:,0] -.2) ).squeeze()
    GMRT_results[z] = n.min(GMRT[z][2:,2])
    print('results: {0},\t{1}'.format(z,n.sqrt(GMRT_results[z])))
    ax.errorbar(float(z), GMRT_results[z], GMRT_results[z]/1.5, fmt='p',ecolor='gray',color='gray', uplims=True, label='Paciga, 2013' if i ==0 else "",capsize = capsize)

#Get MWA 32 data
MWA_results = {}
MWA = eor_results.MWA_32T_all()
print('Results: Z,\t Upper Limits')
print('MWA 32')
for i,z in enumerate(MWA.keys()):
    #index = n.argwhere(MWA[z][:,0] - .2 < .01).squeeze()
    freq= pspec.z2f(z)
    k_horizon = n.sqrt(cosmo_units.eta2kparr(30./c,z)**2 + \
                    cosmo_units.u2kperp(15*freq*1e6/c,z)**2)
    index = n.argwhere( abs(MWA[z][:,0] > k_horizon) ).squeeze()
    MWA_results[z] = n.min(MWA[z][index,2])
    print('results: {0},\t{1}'.format(z,n.sqrt(MWA_results[z])))
    ax.errorbar(float(z), MWA_results[z], MWA_results[z]/1.5, fmt = 'r*', uplims=True, label='Dillon, 2014' if i ==0 else "",capsize = capsize)

MWA128_results = {}
MWA128 = eor_results.MWA_128_all()
print('MWA 128 Beardsley')
for i,z in enumerate(MWA128.keys()):
    #index = n.argwhere(MWA[z][:,0] - .2 < .01).squeeze()
    index = n.argmin(abs(MWA128[z][:,0] - .2)).squeeze()
    MWA128_results[z] = n.min(MWA128[z][index,2])
    print('results: {0},\t{1}'.format(z,n.sqrt(MWA128_results[z])))
    ax.errorbar(float(z), MWA128_results[z], MWA128_results[z]/1.5, fmt = 'g*', uplims=True, label='Dillon, 2015' if i ==0 else "",capsize = capsize)

MWA128_results = {}
MWA128 = eor_results.MWA_128_beardsley_2016_all(pol='NS')
print('MWA 128')
for i,z in enumerate(MWA128.keys()):
    #index = n.argwhere(MWA[z][:,0] - .2 < .01).squeeze()
    index = n.argmin(abs(MWA128[z][:,0] - .2)).squeeze()
    MWA128_results[z] = n.min(MWA128[z][index,2])
    print('results: {0},\t{1}'.format(z,n.sqrt(MWA128_results[z])))
    ax.errorbar(float(z), MWA128_results[z], MWA128_results[z]/1.5, fmt = 'b.', uplims=True, label='Beardsley, 2016' if i ==0 else "",capsize = capsize)
print('LOFAR Patil')
LOFAR_Patil = eor_results.LOFAR_Patil_2017()
print LOFAR_Patil[8.3].shape
Patil_z,Patil_pspec = eor_results.k_slice(0.12,LOFAR_Patil)
print('results: z=',Patil_z,'pspec = ',Patil_pspec)
ax.errorbar(Patil_z,Patil_pspec[:,2],Patil_pspec[:,2]/1.5,fmt='bd',uplims=True,label='Patil, 2017',capsize=capsize)

#Get Paper-32 data
PSA32 = eor_results.PAPER_32_all()
PSA32_results = {}
Jacobs_et_al=[0,1,2,4]
print('PSA32')
for i,z in enumerate( PSA32.keys()):
    index = n.argmin( abs(PSA32[z][:,0] - .2) ).squeeze()
    PSA32_results[z] = n.min(PSA32[z][index,2])
    print('results: {0},\t{1}'.format(z,n.sqrt(PSA32_results[z])))

    if i in Jacobs_et_al and plot_psa32_multiz:
        ax.errorbar(float(z), PSA32_results[z], PSA32_results[z]/1.5, fmt = 'md', uplims=True, label='Jacobs, 2015' if i ==0 else "",capsize = capsize)
    elif plot_psa32 and i not in Jacobs_et_al:
        ax.errorbar(float(z), PSA32_results[z], PSA32_results[z]/1.5, fmt = 'cv', uplims=True, label='Parsons, 2014',capsize = capsize)

#Get PAPER-64 results
if plot_psa64:
    PSA64 = eor_results.PAPER_64_all()
    PSA64_results = {}
    print('PSA64')
    for z in PSA64.keys():
        index = n.argmin( abs(PSA64[z][:,0] - .2)).squeeze()
        PSA64_results[z] = n.min(abs(PSA64[z][index,2]))
        print('results: {0},\t{1}'.format(z,n.sqrt(PSA64_results[z])))
        ax.errorbar(float(z), PSA64_results[z], PSA64_results[z]/1.5, fmt='bs', uplims=True, label= 'Ali, 2015',capsize = capsize)

#zs = [10.87,8.37]
results= {}
if plot_psa64_multiz:
    print('PSA64 multiz')
    for i,fi in enumerate(files):
        z = zs[i]
        f= n.load(fi)
        results_array = n.array([f['k'],f['k3pk'],f['k3pk']+f['k3err'],f['k3pk']-f['k3err']]).T
        negs = n.argwhere(f['k3pk'] < 0).squeeze()

        try: len(negs)
        except: negs = n.array([negs.item()])
        if len(negs) > 0 :
            results_array[negs] = n.array([f['k'][negs],abs(f['k3pk'][negs]),abs(f['k3pk'][negs])+abs(f['k3err'][negs]),abs(f['k3pk'][negs])-abs(f['k3err'][negs])]).T
        index = n.argmin( abs(results_array[:,0] - .2) ).squeeze()
        results[z] = n.min(abs(results_array[index,2]))
        print('results: {0},\t{1}'.format(z,n.sqrt(results[z])))
        ax.errorbar(float(z), results[z], results[z]/1.5, fmt= 'ko', uplims=True,
            label ='Kolopanis, 2016' if i ==0 else "", capsize = capsize)

ax.set_yscale('log')
ax.set_ylabel('$\Delta^{2} (mK)^{2}$')
ax.set_ylim([1e-1,1e7])
ax.set_xlabel('z')
ax.grid(axis='y')

#Add model data.
if not opts.models is None:
    print 'Plotting 21cmFAST Model'
    simk = 0.2
    xlim = ax.get_xlim() #save the data xlimits
    from twentyonecmfast_tools import load_andre_models,all_and
    from glob import glob
    parm_array,k_array,delta2_array,delta2_err_array = load_andre_models(opts.models)
    k_index = n.abs(k_array[0]-simk).argmin()
    alphaXs = n.sort(list(set(parm_array[:,3])))
    Mmins = n.sort(list(set(parm_array[:,4])))
    Nxs = n.sort(list(set(parm_array[:,2])))
    for Nx in Nxs:
        for alphaX in alphaXs:
            for Mmin in Mmins:
                _slice = n.argwhere(all_and([
                                    parm_array[:,2]==Nx,
                                    parm_array[:,3]==alphaX,
                                    parm_array[:,4]==Mmin]
                                    ))
                if len(_slice)==0:continue
                if alphaX==0:
                    label='Cold Reionization'
                    ls = ':k'
                else:
                    label = "Fiducial 21cmFAST model"
                    ls = '-k'
                ax.plot(parm_array[_slice,0],delta2_array[_slice,k_index],
                    ls)#,label=label)
    ax.set_xlim(xlim) #reset to the data xlimits

handles, labels = ax.get_legend_handles_labels()
print [(h,cnt) for cnt,h in enumerate(handles)]
handles = [h[0] if cnt > (len(parm_array)-1) else h for cnt,h in enumerate(handles)]
#handles.insert(-1,handles.pop(0))
#labels.insert(-1,labels.pop(0))
box = ax.get_position()
ax.set_position([box.x0, box.height * .2 + box.y0, box.width, box.height*.8])
#fig.subplots_adjust(bottom=.275,top=.8)
ax.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.4),ncol=4,fontsize=fontsize/2,markerscale=1)
#ax.legend(loc='bottom',ncol=3)
ax.text(10.3,25,"Fiducial 21cmFAST model",fontsize=textsize)
ax.text(10.4,1.5e3,"Cold Reionization",fontsize=textsize)

#add a HERA sensitivity curve
F = n.load('/Users/djacobs/src/21cmSense/hera127drift_pess_0.150.npz')
freqs = n.linspace(80,200)
zs = 1421./freqs - 1
ax.plot(zs,F['T_errs'][3] * (freqs/180.)**-2.55,color='orange',label='HERA127')
ax.text(10.3,2,"HERA127",color='orange',fontsize=textsize)
fig.savefig('eor_results.png')

if opts.plot:
    p.show()
