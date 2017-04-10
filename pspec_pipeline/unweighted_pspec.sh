#!/bin/bash

#Use this script to run an unweighted psa64 or 128 power spectrum
#usage:
#siglos_bash.sh <input lstbinned data> <output power spectrum folder>

indir=$1
outdir=$2
POL='I'
SEP='sep0,1'
chans='95_115'
even_lsts='lst*242.[3456]*'
odd_lsts='lst*243.[3456]*'
appelation='.uvGAL'

noise=''
boot=60
t_eff=69
bl_length=30
window='none'

declare -A rmbls
rmbls['sep0,1']='15_16,0_26,0_44,16_62,3_10,3_25'

cal=psa6240_v003
scriptsdir=/home/mkolopanis/src/capo

# Building up the output directory file structure
# name_of_run/polarization/channel_range/separation/injection_level

for chan in $chans; do

    chandir=$outdir/$chan
    test -e $chandir || mkdir -p $chandir

    for pol in $POL; do

        poldir=$chandir/$pol
        test -e $poldir || mkdir -p $poldir

        for sep in $SEP; do

            sepdir=$poldir/$sep
            test -e $sepdir || mkdir -p $sepdir

            EVEN_FILES=${indir}'/even/'${sep}/${even_lsts}$appelation
            ODD_FILES=${indir}'/odd/'${sep}/${odd_lsts}$appelation

            ${scriptsdir}/pspec_pipeline/pspec_cov_v004.py --window=${window}\
             -a cross -p $pol -c ${chan} -C ${cal} -b ${boot}\
              ${noise} --rmbls=${rmbls[$sep]}\
               --output=${sepdir} ${EVEN_FILES} ${ODD_FILES}

            echo "${scriptsdir}/pspec_pipeline/sigloss_sim.py\
             --window=${window} -a cross -p ${pol} -c ${chan}\
              -C ${cal} -b ${boot} ${noise} --rmbls=${rmbls}\
               --output=${sepdir} ${EVEN_FILES} ${ODD_FILES} " > ${sepdir}/notes.tx
        done
    done
done

#Run through each combination of power spectra and compute
#signal loss corrected spectrum
for chan in $chans; do

    chandir=$outdir/$chan

    for pol in $POL; do

        poldir=$chandir/$pol

        for sep in $SEP;do

            sepdir=$poldir/$sep

            python ${scriptsdir}/pspec_pipeline/boots_to_unweighted.py\
             --t_eff=${t_eff} --bl_length=${bl_length} --outfile=$sepdir $sepdir

        done
    done
done


for chan in $chans; do

    chandir=$outdir/$chan

    for pol in $POL; do

        poldir=$chandir/$pol

        for sep in $SEP; do

            sepdir=$poldir/$sep

            ${scriptsdir}/mjk/scripts/plot_upper_lims_simple.py\
             $sepdir/pspec_pI_pk_k3pk.npz\
              --noisefiles='/home/mkolopanis/psa64/21cmsense_noise/psa64_09Jan2017/*drift_pess*152.npz'\
               --outfile="${outdir}/pspec_unweighted_${outdir}_${chan}_${pol}_${sep}" #--psa32 --psa32_noise='/home/mkolopanis/psa64/21cmsense_noise/psa32_noise/*drift_mod*1[5]0.npz'
        done
    done
done
