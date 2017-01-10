#!/bin/bash

#Use this script to run a psa64 or 128 power spectrum
#with signal loss corrections
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
inject_range='numpy.logspace(-2,3,10)'

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


            for inject in `python -c "import numpy; print ' '.join(map(str, ${inject_range}))"` ; do

                injdir=$sepdir/"inject_${inject}"
                test -e ${injdir} || mkdir -p ${injdir}
                echo SIGNAL_LEVEL=${inject}


                ${scriptsdir}/pspec_pipeline/sigloss_sim.py --window=${window} -a cross -p $pol -c ${chan} -C ${cal} -b ${boot} -i ${inject} ${noise} --rmbls=${rmbls[$sep]} --output=${injdir} ${EVEN_FILES} ${ODD_FILES}

                echo "${scriptsdir}/pspec_pipeline/sigloss_sim.py --window=${window} -a cross -p ${pol} -c ${chan} -C ${cal} -b ${boot} -i ${inject} ${noise} --rmbls=${rmbls} --output=${injdir} ${EVEN_FILES} ${ODD_FILES} " > ${injdir}/notes.txt
            done
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

            python ${scriptsdir}/pspec_pipeline/boots_to_pspec.py --t_eff=${t_eff} --bl_length=${bl_length} --outfile=$sepdir $sepdir
            python /${scriptsdir}/pspec_pipeline/sigloss_limits.py --outfile=$sepdir $sepdir/inject_*/pspec_pk_k3pk*.npz

        done
    done
done


for chan in $chans; do

    chandir=$outdir/$chan

    for pol in $POL; do

        poldir=$chandir/$pol

        for sep in $SEP; do

            sepdir=$poldir/$sep

            ${scriptsdir}/mjk/scripts/plot_upper_lims_simple.py    $sepdir/pspec_limits_k3pk_p[CI]_85.npz --noisefiles='/home/mkolopanis/psa64/21cmsense_noise/dish_size_1/*drift_mod*150.npz'   --outfile="${outdir}/pspec_${outdir}_${chan}_${pol}_${sep}" #--psa32 --psa32_noise='/home/mkolopanis/psa64/21cmsense_noise/psa32_noise/*drift_mod*1[5]0.npz'
        done
    done
done
