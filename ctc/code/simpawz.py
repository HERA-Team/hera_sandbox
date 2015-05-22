#! /usr/bin/env python

"""
NAME: 
        simpawz.py
AUTHOR: 
        Carina Cheng
PURPOSE: 
        This code can create maps of foregrounds using the Global Sky Model and/or maps of a theoretical cosmological signal, and can simulate visibilities of the maps
        Run this code on Folio!
SCRIPTS NEEDED:
        gsmmf.sh
        gsmtofits.py
        pspec_sim_v3.py
        batch_vissim.sh
        vis_simulation_v4.py
        combine_times.py
        vis_simulation_K2Jy.py
OUTPUTS:
        GSM: maps labeled gsm1001.fits, gsm1002.fits, etc.
        PSPEC: maps labeled pspec1001.fits, pspec1002.fits, etc.
        SIM: concatenated UV files labeled pspec_Jy.uv/gsm_Jy.uv [Jy] and pspec_K.uv/gsm_K.py [K]
EXAMPLE CALL:
        ./simpawz.py --gsm --pspec --sim
"""


import os, sys
import aipy
import numpy
import scipy
import glob
import optparse
import subprocess
import time
import datetime

#GLOBAL VARIABLES                                ### to be updated each time ###

NCHAN = 5 #number of frequency channels
SFREQ = 0.1 #starting frequency [GHz]
SDF = 0.001 #spacing between frequencies [GHz]

INTTIME = 20000 #integration time [s]
STARTJD = 2454500 #starting julian date
ENDJD = 2454501 #ending julian date

FILEPATH = '/home/cacheng/capo/ctc/simpawz_test/' #path to save all outputs (either where your maps are or where you want them to be created)

LMAX = 100 #maximum l to simulate PSPEC maps up to (if opts.pspec == True)
K_VALS = numpy.arange(0.001,0.5,0.01) #k-values to make maps for (if opts.pspec == True)
PK_VALS = 0.000505*(2*numpy.pi**2)/(K_VALS**3) #Pk-values to make maps for (if opts.pspec == True)

ANT = '0_16' #antenna numbers for baseline simulated (if opts.vis == True)
CAL = 'psa898_v003' #calfile path (if opts.vis == True)
NJOBS = 3 #number of jobs on Folio to run (if opts.vis == True)


#OPTIONS

o = optparse.OptionParser()
o.set_usage('simpawz.py [options]')
o.set_description(__doc__)
o.add_option('--gsm', action='store_true', default=False, help='Make GSM maps.')
o.add_option('--pspec', action='store_true', default=False, help='Make PSPEC maps.')
o.add_option('--vis', action='store_true', default=False, help='Simulate visibilities.')
opts, args = o.parse_args(sys.argv[1:])

os.chdir(FILEPATH)


#GSM

if opts.gsm == True:

    proc = subprocess.Popen(["which gsmmf.sh"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    os.chdir(out[:-9])
    gsmargs = open('args.dat','w')
    gsmargs.write(' '.join(['gsm',str(SFREQ*1000),str(SDF*1000),str(NCHAN)]))
    gsmargs.close()
    os.chdir(FILEPATH)
       
    if os.path.exists('gsm1001.fits') == True:
        if opts.vis == False:
            ans = raw_input('GSM .fits files already exist. Delete and re-make the maps [y/n]? ')
        if opts.vis == True:
            ans = raw_input('GSM .fits files already exist. Delete and re-make the maps [y] or go straight to simulating visibilities [n]? ')
        if ans == 'y':
            print 'Deleting GSM .fits files...'
            os.system("rm -r gsm*fits")
            print 'Making GSM maps...'
            os.chdir(out[:-9])
            os.system(out) #makes maps
            os.system('mv gsm*dat '+FILEPATH)
            os.chdir(FILEPATH)

            for name in glob.glob(os.path.join('*.dat')):
                d = numpy.loadtxt(name)
                h = aipy.healpix.HealpixMap(nside=512)
                h.map = d
                h.to_fits(name.replace('.dat','.fits')) #converts to fits files
                print name + ' -> ' + name.replace('.dat','.fits')
        print 'GSM maps have been made!'
 
    if os.path.exists('gsm1001.dat') == False:

        print 'Making GSM maps...'
        os.chdir(out[:-9])
        os.system(out) #makes maps
        os.system('mv gsm*dat '+FILEPATH) #moves sky mapsi
        os.chdir(FILEPATH)

        for name in glob.glob(os.path.join('*.dat')):
            d = numpy.loadtxt(name)
            h = aipy.healpix.HealpixMap(nside=512)
            h.map = d
            h.to_fits(name.replace('.dat','.fits')) #converts to fits files
            print name + ' -> ' + name.replace('.dat','.fits')
            
        print 'GSM maps have been made!'


#PSPEC

if opts.pspec == True:

    k_vals_string = ' '.join(map(str,K_VALS))
    Pk_vals_string = ' '.join(map(str,PK_VALS))
    k_and_Pk = k_vals_string+' '+Pk_vals_string
    
    if os.path.exists('pspec1001.fits') == True:
        if opts.vis == False:
            ans = raw_input('PSPEC .fits files already exist. Delete and re-make the maps [y/n]? ')
        if opts.vis == True:
            ans = raw_input('PSPEC .fits files already exist. Delete and re-make the maps [y] or go straight to simulating visibilities [n]? ')
        if ans == 'y':
            print 'Deleting PSPEC .fits files...'
            os.system("rm -r pspec*fits")
            print 'Making PSPEC maps...'
            command = 'pspec_sim_v3.py --nchan '+str(NCHAN)+' --sdf '+str(SDF)+' --sfreq '+str(SFREQ)+' --lmax '+str(LMAX)+' --path '+FILEPATH+' '+k_and_Pk
            os.system(command)
            print 'PSPEC maps have been made!'
 
    if os.path.exists('pspec1001.fits') == False:

        print 'Making PSPEC maps...'
        command = 'pspec_sim_v3.py --nchan '+str(NCHAN)+' --sdf '+str(SDF)+' --sfreq '+str(SFREQ)+' --lmax '+str(LMAX)+' --path '+FILEPATH+' '+k_and_Pk
        os.system(command)
        print 'PSPEC maps have been made!'

       
#VIS 

if opts.vis == True:

    if opts.gsm == True:

        print 'Simulation GSM visibility...'
 
        shellscript = open('batch_sim.sh','w')
        shellscript.write('#$ -S /bin/bash'+'\n'+'#$ -V'+'\n'+'#$ -cwd'+'\n'+'#$ -l h_vmem=16G'+'\n'+'#$ -l paper'+'\n'+'#$ -o /data2/home/cacheng/capo/ctc/code/gridoutput'+'\n'+'#$ -e /data2/home/cacheng/capo/ctc/code/gridoutput')
        shellscript.write('\n\n')
        shellscript.write('myargs=`pull_args.py $*`'+'\n\n')
        shellscript.write('echo my times: ${myargs}'+'\n\n')
        shellscript.write('name=`echo ${myargs} | cut -d " " -f 1`'+'\n')
        shellscript.write('echo first arg: ${name}'+'\n\n')

        command = 'vis_simulation_v4.py --sdf '+str(SDF)+' --sfreq '+str(SFREQ)+' --nchan '+str(NCHAN)+' --inttime '+str(INTTIME)+' --map gsm --mappath '+FILEPATH+' --filename '+FILEPATH+'gsm_${name}.uv -C '+CAL+' -a '+ANT+' ${myargs}'

        shellscript.write('echo '+command+'\n\n')
        shellscript.write(command)
        shellscript.close()

        qsubcommand = 'qsub -t 1-'+str(NJOBS)+' batch_sim.sh `python -c "import numpy; import aipy; print ' + """'"""+""" '"""+'.join(map(str,numpy.arange('+str(STARTJD)+','+str(ENDJD)+','+str(INTTIME)+'/aipy.const.s_per_day)))"`'
       
        os.system(qsubcommand) #runs simulation code
 
        running=True
        while running==True:
            print 'Checking simulation status at '+str(datetime.datetime.now())
            status = subprocess.Popen('qstat | grep cacheng',stdout=subprocess.PIPE,shell=True) 
            (out,err) = status.communicate()
            out = str(out)
            if out.find("cacheng") != -1:
                print '   Simulation still running. Checking again in another 30 seconds.'
                running=True
                time.sleep(30) #checks every 30 seconds
            else:
                running=False

        print 'UV files made!'
        print 'Combining times now...'
    
        command = 'combine_times.py --uvnew gsm_K.uv gsm*.uv'
        os.system(command)

        print 'Converting from [K] to [Jy]...'

        command = 'vis_simulation_K2Jy.py --uvold gsm_K.uv --nchan '+str(NCHAN)+' --sdf '+str(SDF)+' --sfreq '+str(SFREQ)+' --uvnew gsm_Jy.uv'
        os.system(command)

        print 'Done! Final UV file is gsm_Jy.uv.'
   
    if opts.pspec == True:
         
        shellscript = open('batch_sim.sh','w')
        shellscript.write('#$ -S /bin/bash'+'\n'+'#$ -V'+'\n'+'#$ -cwd'+'\n'+'#$ -l h_vmem=16G'+'\n'+'#$ -l paper'+'\n'+'#$ -o /data2/home/cacheng/capo/ctc/code/gridoutput'+'\n'+'#$ -e /data2/home/cacheng/capo/ctc/code/gridoutput')
        shellscript.write('\n\n')
        shellscript.write('myargs=`pull_args.py $*`'+'\n\n')
        shellscript.write('echo my times: ${myargs}'+'\n\n')
        shellscript.write('name=`echo ${myargs} | cut -d " " -f 1`'+'\n')
        shellscript.write('echo first arg: ${name}'+'\n\n')
       
        command = 'vis_simulation_v4.py --sdf '+str(SDF)+' --sfreq '+str(SFREQ)+' --nchan '+str(NCHAN)+' --inttime '+str(INTTIME)+' --map pspec --mappath '+FILEPATH+' --filename '+FILEPATH+'pspec_${name}.uv -C '+CAL+' -a '+ANT+' ${myargs}'

        shellscript.write('echo '+command+'\n\n')
        shellscript.write(command)
        shellscript.close()
  
        qsubcommand = 'qsub -t 1-'+str(NJOBS)+' batch_sim.sh `python -c "import numpy; import aipy; print ' + """'"""+""" '"""+'.join(map(str,numpy.arange('+str(STARTJD)+','+str(ENDJD)+','+str(INTTIME)+'/aipy.const.s_per_day)))"`'
      
        os.system(qsubcommand)
        
        running=True
        while running==True:
            print 'Checking simulation status at '+str(datetime.datetime.now())
            status = subprocess.Popen('qstat | grep cacheng',stdout=subprocess.PIPE,shell=True) 
            (out,err) = status.communicate()
            out = str(out)
            if out.find("cacheng") != -1:
                print '   Simulation still running. Checking again in another 30 seconds.'
                running=True
                time.sleep(30) #checks every 30 seconds
            else:
                running=False

        print 'UV files made!'
        print 'Combining times now...'
    
        command = 'combine_times.py --uvnew pspec_K.uv pspec*.uv'
        os.system(command)

        print 'Converting from [K] to [Jy]...'
    
        command = 'vis_simulation_K2Jy.py --uvold pspec_K.uv --nchan '+str(NCHAN)+' --sdf '+str(SDF)+' --sfreq '+str(SFREQ)+' --uvnew pspec_Jy.uv'
        os.system(command)

        print 'Done! Final UV file is pspec_Jy.uv.'
        
