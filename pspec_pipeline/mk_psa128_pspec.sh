#! /bin/bash
export PYTHONPATH='.':$PYTHONPATH

(
echo using config $*
. $*
#set defaults to parameters which might not be set
if [[ ! -n ${WINDOW} ]]; then export WINDOW="none"; fi
#for posterity, print the cal file we are using

pywhich $cal

threadcount=`python -c "c=map(len,['${pols}'.split(),'${chans}'.split(),'${seps}'.split()]);print c[0]*c[1]*c[2]"`
echo Running $threadcount pspecs
PIDS=""

#FILES=`lst_select.py -C ${cal} --ra=${RA} ${DATAPATH}`
test -e ${PREFIX} || mkdir $PREFIX
for chan in $chans; do
    chandir=${PREFIX}/${chan}
    test -e ${chandir} || mkdir ${chandir}
    for pol in $pols; do
        echo "Starting work on ${pol}" 
        poldir=${chandir}/${pol}
        test -e ${poldir} || mkdir ${poldir}
        if [ ! -e ${poldir}/pspec_${PREFIX}_${chan}_${pol}.png ]; then
            for sep in $seps; do
                sepdir=${poldir}/${sep}
                #form up the path to the data use ()s for globing
                EVEN_FILES=(${EVEN_DATAPATH}/${sep}/*${FILEAPPELLATION})
                #convert from an array to a... list? ida know. bash stuff.
                EVEN_FILES=`lst_select.py -C ${cal} --ra=${LST} ${EVEN_FILES[@]}`
                ODD_FILES=(${ODD_DATAPATH}/${sep}/*${FILEAPPELLATION})
                ODD_FILES=`lst_select.py -C ${cal} --ra=${LST} ${ODD_FILES[@]}`
                test -e ${sepdir} || mkdir ${sepdir}
                LOGFILE=`pwd`/${PREFIX}/${chan}_${pol}_${sep}.log
                echo this is mk_pspec.sh with  |tee -a  ${LOGFILE}
                echo recording to ${LOGFILE} | tee -a ${LOGFILE}
                echo experiment: ${PREFIX}|tee -a ${LOGFILE}
                echo channels: ${chan}|tee -a ${LOGFILE}
                echo polarization: ${pol}|tee -a ${LOGFILE}
                echo separation: ${sep}|tee -a ${LOGFILE}
                echo `date` | tee -a ${LOGFILE}

                #ANTS=`grid2ant.py -C ${cal} --seps="${sep}"`
                ANTS='cross'
                echo Beginning pspec calculation | tee -a ${LOGFILE}
                echo using ${#EVEN_FILES} even files and ${#ODD_FILES} odd files
                echo python ${SCRIPTSDIR}/pspec_cov_v002.py -C ${cal} \
                     -b ${NBOOT} -a ${ANTS} -c ${chan} -p ${pol}\
                      --window=${WINDOW}  ${NOPROJ} --output=${sepdir} \
                       ${EVEN_FILES} ${ODD_FILES} ${OPTIONS}
                
                python ${SCRIPTSDIR}/pspec_cov_v003.py -C ${cal} -b ${NBOOT} \
                    -a ${ANTS} -c ${chan} -p ${pol} --window=${WINDOW} \
                      ${NOPROJ} --output=${sepdir} \
                      ${EVEN_FILES} ${ODD_FILES} ${OPTIONS} #\
                     #| tee -a ${LOGFILE}
                if [ $? -ne 0 ] 
                then
                exit
                fi
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
        mv  pspec_${PREFIX}_${chan}_${pol}.png ${poldir}/
        mv  posterior_${PREFIX}_${chan}_${pol}.png ${poldir}/
        mv posterior.txt ${poldir}/
        mv pspec_${PREFIX}_${chan}_${pol}.npz ${poldir}/
    done
done
)
