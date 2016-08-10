#! /bin/bash

export PYTHONPATH='.':$PYTHONPATH
(
echo using config $*
. $*
#set defaults to parameters which might not be set
if [[ ! -n ${WINDOW} ]]; then export WINDOW="none"; fi
#for posterity Print the cal file we are using

pywhich $cal

threadcount=`python -c "c=map(len,['${pols}'.split(),'${chans}'.split(),'${seps}'.split()]);print c[0]*c[1]*c[2]"`
echo Running $threadcount pspecs
PIDS=""



test -e ${PREFIX} || mkdir $PREFIX
for chan in $chans; do
    chandir=${PREFIX}/${chan}
    test -e $chandir || mkdir $chandir
    for pol in ${pols}; do
        echo "Starting ${pol}"
        poldir=${chandir}/${pol}
        test -e ${poldir} || mkdir ${poldir}
        if [ ! -e ${poldir}/pspec_${PREFIX}_${chan}_${pol}.png ]; then
            for sep in $seps; do
               sepdir=${poldir}/${sep}
               EVEN_FILES=${EVEN_DATAPATH}${sep}/*242.[3456]*uvGL
               ODD_FILES=${ODD_DATAPATH}${sep}/*243.[3456]*uvGL
               test -e ${sepdir} || mkdir ${sepdir}
               LOGFILE=`pwd`/${PREFIX}/${chan}_${pol}_${sep}.log
               echo this is make_psa64_pspec_cav.sh with  |tee  ${LOGFILE}
               echo experiment: ${PREFIX}|tee -a ${LOGFILE}
               echo channels: ${chan}|tee -a ${LOGFILE}
               echo polarization: ${pol}|tee -a ${LOGFILE}
               echo separation: ${sep}|tee -a ${LOGFILE}
               echo `date` | tee -a ${LOGFILE} 
               ANTS='cross'
                   
               echo python ${SCRIPTSDIR}/pspec_cov_full_band.py \
                    --window=${WINDOW} -p ${pol} -c ${chan} --band=${band}\
                     -b ${NBOOT} -C ${cal} --auto\
                     -a ${ANTS} --output=${sepdir} \
                     ${EVEN_FILES} ${ODD_FILES}

               python ${SCRIPTSDIR}/pspec_cov_full_band.py \
                    --window=${WINDOW} -p ${pol} -c ${chan} --band=${band}\
                     -b ${NBOOT} -C ${cal} --auto\
                     -a ${ANTS} --output=${sepdir} \
                     ${EVEN_FILES} ${ODD_FILES} \
                     |tee -a ${LOGFILE}
                
                
                echo beginning bootstrap: `date` | tee -a ${LOGFILE} 
                ${SCRIPTSDIR}/pspec_cov_boot.py ${sepdir}/pspec_boot*npz | tee -a ${LOGFILE} 
                echo complete! `date`| tee -a ${LOGFILE} 
                mv pspec.npz ${sepdir}/
                PIDS="${PIDS} "$!
            done
        fi
    done
done

echo waiting on `python -c "print len('${PIDS}'.split())"` power spectra threads ${PIDS} 
wait $PIDS
echo power spectrum complete

echo averaging power spectra for pols/channels
for chan in $chans; do
    chandir=${PREFIX}/${chan}
    for pol in $pols; do
        echo "Generating plots for ${chan}: ${pol}"
        poldir=${chandir}/${pol}
        #PLOT
        ${SCRIPTSDIR}/plot_pk_k3pk_zsa_2.py ${poldir}/*/pspec.npz 
        mv pspec_pk_k3pk.npz pspec_${PREFIX}_${chan}_${pol}.npz
        mv pspec.png pspec_${PREFIX}_${chan}_${pol}.png
        mv posterior.png posterior_${PREFIX}_${chan}_${pol}.png
        cp  pspec_${PREFIX}_${chan}_${pol}.png ${poldir}/
        cp  posterior_${PREFIX}_${chan}_${pol}.png ${poldir}/
        cp pspec_${PREFIX}_${chan}_${pol}.npz ${poldir}/
    done
done    

    )
