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


#GLOBAL VARIABLES

NCHAN = 5 #number of frequency channels
SFREQ = 0.1 #starting frequency [GHz]
SDF = 0.001 #spacing between frequencies [GHz]

INTTIME = 20000 #integration time [s]
STARTJD = 2454500 #starting julian date
ENDJD = 2454501 #ending julian date

CAL = 'psa898_v003' #calfile path

LMAX = 100 #maximum l to simulate PSPEC maps up to
K_VALS = numpy.arange(0.001,0.5,0.01) #k-values to make maps for
PK_VALS = 0.000505*(2*numpy.pi**2)/(K_VALS**3) #Pk-values to make maps for

ANT = '0_16'

FILEPATH = '/home/cacheng/capo/ctc/simpawz_test/' #path to save all outputs


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

    print 'Simulating visibility...'

    shellscript = open('batch_vissim2.sh','w')
    shellscript.write('#$ -S /bin/bash'+'\n'+'#$ -V'+'\n'+'#$ -cwd'+'\n'+'#$ -l h_vmem=16G'+'\n'+'#$ -l paper'+'\n'+'#$ -o /data2/home/cacheng/capo/ctc/code/gridoutput'+'\n'+'#$ -e /data2/home/cacheng/capo/ctc/code/gridoutput')
    shellscript.write('\n\n')
    shellscript.write('myargs=`pull_args.py $*`'+'\n\n')
    shellscript.write('echo my times: ${myargs}'+'\n\n')
    shellscript.write('name=`echo ${myargs} | cut -d " " -f 1`'+'\n')
    shellscript.write('echo first arg: ${name}'+'\n\n')

    if opts.gsm == True:
        
        command = 'vis_simulation_v4.py --sdf '+str(SDF)+' --sfreq '+str(SFREQ)+' --nchan '+str(NCHAN)+' --inttime '+str(INTTIME)+' --map gsm --mappath '+FILEPATH+' --filename '+FILEPATH+'gsm_${name}.uv -C '+CAL+' -a '+ANT+' ${myargs}'

        shellscript.write('echo '+command+'\n\n')
        shellscript.write(command)
        shellscript.close()
        
    if opts.pspec == True:
        
        command = 'vis_simulation_v4.py --sdf '+str(SDF)+' --sfreq '+str(SFREQ)+' --nchan '+str(NCHAN)+' --inttime '+str(INTTIME)+' --map pspec --mappath '+FILEPATH+' --filename '+FILEPATH+'pspec_${name}.uv -C '+CAL+' -a '+ANT+' ${myargs}'

        shellscript.write('echo '+command+'\n\n')
        shellscript.write(command)
        shellscript.close()
  

    
