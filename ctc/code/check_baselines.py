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

o = optparse.OptionParser()
o.set_usage('check_baselines.py [options] *npz')
o.set_description(__doc__)
o.add_option('-a',dest='ants',default='all',
            help='List of baselines (ex: 64_49,3_25). Default is "all".')
o.add_option('--plot',dest='plot',default=False,action="store_true",
            help='Plot LST model for each baseline.')
o.add_option('--good',dest='good',default='64_49',
            help='Good baseline to use as comparison. Default is 64_49.')
opts,args = o.parse_args(sys.argv[1:])


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

lendata = len(firstfile[bls[0]])

gooda1,gooda2 = opts.good.split('_')[0],opts.good.split('_')[1]
bl_good = str(aipy.miriad.ij2bl(int(gooda1),int(gooda2)))
bl_ant_good = '('+str(gooda1)+','+str(gooda2)+')'

good_ind = numpy.where(numpy.array(bls) == bl_good)[0][0]
bls_ants.insert(0,bls_ants.pop(good_ind)) #brings good baseline to front of baseline array
bls.insert(0,bls.pop(good_ind))

### All baselines in the npz slices are 30-m E/W baselines ###
### There are 98 of them ###

for bl in range(len(bls)): #loop over baselines

    mybl = bls[bl]
    print str(bl+1)+'/'+str(len(bls))+': reading baseline '+str(bls_ants[bl])

    all_data = []
    all_lsts = []

    for filename in args: #loop over files to get data and LSTs
        jd_int = int(filename.split('.')[1])
        jd_dec = float(filename.split('.')[1]+'.'+filename.split('.')[2])
        npz = numpy.load(filename)
        try:
            data = npz[mybl]
            lsts = npz['t'+mybl]*12/numpy.pi #LST hours  
            all_data.append(data)
            all_lsts.append(lsts)
        except:
            continue

    all_data = numpy.array(all_data)
    all_lsts = numpy.array(all_lsts)
    all_data = numpy.concatenate(all_data)
    all_lsts = numpy.concatenate(all_lsts)

    #bin LSTs and plot per baseline
    all_lsts,all_data = zip(*sorted(zip(all_lsts,all_data)))
    num_bins = 100
    d_lst = float((all_lsts[-1]-all_lsts[0])/num_bins)
    lst_grid = numpy.arange(all_lsts[0],all_lsts[-1]+d_lst,d_lst)
    model = numpy.zeros(num_bins)
    avg_lsts = numpy.zeros(num_bins)
    for i in range(num_bins):
        bin = []
        for j in range(len(all_lsts)):
            if all_lsts[j] >= lst_grid[i] and all_lsts[j] < lst_grid[i+1]:
                bin.append(all_data[j])
        med_val = numpy.median(bin)
        model[i] = med_val.real
        avg_lsts[i] = (lst_grid[i+1]+lst_grid[i])/2
    if bl == 0:
        compare_model = model

    #If model deviates from compare_model significantly, print out the name of the baseline here

    #plotting LST model
    if opts.plot == True:
        plt.subplot(121)
        plt.plot(avg_lsts,model,label=bls_ants[bl])
        plt.xlabel('LST')
        plt.ylabel('Jy')
        plt.subplot(122)
        plt.plot(avg_lsts,model/compare_model,label=bls_ants[bl])    
        plt.xlabel('LST')
        plt.ylabel('Normalized to '+str(bl_ant_good))

if opts.plot == True:   
    plt.legend(loc=1,prop={'size':8})
    plt.suptitle('LST Models for Different Baselines')
    plt.show()



