#! /usr/bin/env python

###
# NAME:    check_baselines.py
# AUTHOR:  Carina Cheng
# PURPOSE: Extracts a certain type of baseline from PSA128 slice data (npz files) and compares LST models between them
###

import numpy
import matplotlib.pyplot as plt
import aipy
import os, sys
import optparse


#options
o = optparse.OptionParser()
o.set_usage('check_baselines.py [options] *npz')
o.set_description(__doc__)
o.add_option('-a',dest='ants',default='all',
            help='List of baselines (ex: 64_49,3_25). Default is "all".')
o.add_option('--models',dest='models',default=False,action="store_true",
            help='Plot LST model for each baseline.')
o.add_option('--grids',dest='grids',default=False,action="store_true",
            help='Plot grid of data for each baseline.')
o.add_option('--res',dest='res',default=False,action="store_true",
            help='Plot (data-LST model) for each baseline.')
opts,args = o.parse_args(sys.argv[1:])


#load baselines
bls = [] #baselines by aipy convention
bls_ants = [] #baselines by antenna numbers
firstfile = numpy.load(args[0])

if opts.ants != 'all':
    antennas = opts.ants.split(',')
    for a in range(len(antennas)):
        a1,a2 = antennas[a].split('_')[0],antennas[a].split('_')[1]
        ant_string = '('+str(a1)+','+str(a2)+')'
        bls_ants.append(ant_string)
        bls.append(str(aipy.miriad.ij2bl(int(a1),int(a2))))
else:
    keys = firstfile.keys()
    for i in range(len(keys)):
        if keys[i][0] != 't':
            bls.append(keys[i])
            bls_ants.append(aipy.miriad.bl2ij(keys[i]))

#bls = bls[:30] #plot first 30 baselines, for testing purposes
#bls_ants = bls_ants[:30]

lendata = len(firstfile[bls[0]])

plot_grid_num = 1
subplot_grid_count = 1
plot_res_num = 1
subplot_res_count = 1

### All baselines in the npz slices are 30-m E/W baselines ###
### There are 98 of them ###


#get data
for bl in range(len(bls)): #loop over baselines
    mybl = bls[bl]
    print str(bl+1)+'/'+str(len(bls))+': reading baseline '+str(bls_ants[bl])
    all_data = []
    all_lsts = []
    all_jds = []
    for filename in args: #loop over files to get data and LSTs
        jd_int = int(filename.split('.')[1])
        jd_dec = float(filename.split('.')[1]+'.'+filename.split('.')[2])
        npz = numpy.load(filename)
        try:
            data = npz[mybl]
            lsts = npz['t'+mybl]*12/numpy.pi #LST hours  
            all_data.append(data)
            all_lsts.append(lsts)
            for j in range(len(lsts)):
                all_jds.append(jd_int)
        except:
            continue
    all_data = numpy.array(all_data)
    all_lsts = numpy.array(all_lsts)
    all_jds = numpy.array(all_jds)
    all_data = numpy.concatenate(all_data)
    all_lsts = numpy.concatenate(all_lsts)
    
    #bin LSTs and JDs and make grid
    all_lsts,all_data,all_jds = zip(*sorted(zip(all_lsts,all_data,all_jds)))
    num_bins = 100
    d_lst = float((all_lsts[-1]-all_lsts[0])/num_bins)
    lst_grid = numpy.arange(all_lsts[0],all_lsts[-1]+d_lst,d_lst)
    jd_grid = numpy.arange(numpy.min(all_jds),numpy.max(all_jds)+1,1)
    lst_inds = numpy.digitize(all_lsts,lst_grid,right=True)
    jd_inds = numpy.digitize(all_jds,jd_grid,right=True)
    gridded_data = numpy.zeros((len(lst_grid),len(jd_grid))).astype(numpy.complex64)
    for i in range(len(all_lsts)):
        gridded_data[lst_inds[i],jd_inds[i]] = all_data[i]
    gridded_data = numpy.ma.masked_where(gridded_data==0,gridded_data) 
   
    #plot grid 
    if opts.grids == True:
        if subplot_grid_count == 26:
            plot_grid_num += 1
            subplot_grid_count = 1
        plt.figure(plot_grid_num,figsize=(10,10))
        plt.subplot(5,5,subplot_grid_count)
        plt.imshow(numpy.abs(gridded_data),aspect='auto',vmax=0.05,interpolation='nearest',extent=(0,jd_grid.max()-jd_grid.min(),lst_grid.max(),lst_grid.min()))
        #plt.colorbar()
        plt.text(0.92,-0.07,"+%i"%jd_grid.min(),fontsize=6,transform=plt.gca().transAxes)
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.xlabel('JD',fontsize=8)
        plt.ylabel('LST',fontsize=8)
        plt.title(bls_ants[bl],fontsize=10)
        plt.tight_layout()
        subplot_grid_count += 1
    
    #make LST model
    model = numpy.ma.median(gridded_data,axis=1)
    model.shape += (1,)
    
    #plot model
    if opts.models == True:
        fignum = (numpy.ceil(len(bls_ants)/25.)*2+1)
        fig2 = plt.figure(fignum)
        plt.plot(lst_grid,numpy.real(model[:,0]),label=bls_ants[bl])
        plt.xlabel('LST')
        plt.ylabel('Jy')
        plt.legend(loc=1,prop={'size':8})
        plt.title('LST Models for Different Baselines')

    #plot residuals
    if opts.res == True:   
        fignum = numpy.ceil(len(bls_ants)/25.)+plot_res_num
        if subplot_res_count == 26:
            plot_res_num += 1 
            fignum = numpy.ceil(len(bls_ants)/25.)+plot_res_num 
            subplot_res_count = 1 
        plt.figure(fignum,figsize=(10,10))
        plt.subplot(5,5,subplot_res_count)
        plt.imshow(numpy.abs(gridded_data-model),aspect='auto',vmax=0.05,interpolation='nearest',extent=(0,jd_grid.max()-jd_grid.min(),lst_grid.max(),lst_grid.min()))
        #plt.colorbar()
        plt.text(0.92,-0.07,"+%i"%jd_grid.min(),fontsize=6,transform=plt.gca().transAxes)
        plt.tick_params(axis='both', which='major', labelsize=6)
        plt.xlabel('JD',fontsize=8)
        plt.ylabel('LST',fontsize=8)
        plt.title(bls_ants[bl],fontsize=10)
        plt.tight_layout()
        subplot_res_count += 1

#show plots
if opts.grids == True:
    plt.show()
if opts.res == True:
    plt.show()
if opts.models == True:   
    plt.show()



