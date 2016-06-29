#! /bin/bash
export PYTHONPATH='.':$PYTHONPATH

(
echo using config $*
. $*
#set defaults to parameters which might not be set

if [[ ! -n ${WINDOW} ]]; then export WINDOW="none"; fi

#for posterity,
echo the cal file we are using: 
pywhich $cal


#form up the sigloss factors
cov_array=($covs)
chan_array=($chans)

declare -A scales
for (( i=0; i<${#chan_array[@]}; ++i )); do
    if [ -z "$covs" ]; then
        scales[${chan_array[$i]}]=1
    else
        scales[${chan_array[$i]}]=${cov_array[$i]}
    fi
done

##create the proper flags in pspec_cov for the desired injected noise
n_array=($noise)
declare -A noise_append
for (( i=0; i<${#chan_array[@]}; ++i )); do

    if [ -z "$noise" ]; then
        noise_append[${chan_array[$i]}]=""
        continue
    else
        noise_append[${chan_array[$i]}]="--noise=${n_array[$i]}"
    fi

    if [ -n "$teff" ]; then
        noise_append[${chan_array[$i]}]+=" --teff=${teff}"
    fi

    if [ $NOISE_ONLY == True ]; then
        noise_append[${chan_array[$i]}]+=" --noise_only"
    fi

    if [ $FILTER_NOISE == True ]; then
        noise_append[${chan_array[$i]}]+=" --filter_noise"
    fi

    if [ $BOXCAR == True ]; then
        noise_append[${chan_array[$i]}]+=" --boxcar"
    fi

    if [ ! -z "$BL_SCALE" ]; then
        noise_append[${chan_array[$i]}]+=" --bl_scale=${BL_SCALE}"
    fi
    if [ ! -z "$FR_WIDTH_SCALE" ]; then
        noise_append[${chan_array[$i]}]+=" --fr_width_scale=${FR_WIDTH_SCALE}"
    fi
done

if [ $USE_pI == True ]; then
    nocov='--nocov'
else
    nocov=''
fi


#for (( i=0; i<${#chan_array[@]}; ++i )); do
#  echo ${chan_array[$i]} ${noise_append[${chan_array[$i]}]}
#done



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

                if [ $COV == True ]; then
                    #check if files are manually globbed before lst_select 
                    if [[ -z $EVEN_GLOB ]]; then
                        #form up the path to the data use ()s for globing
                        EVEN_FILES=(${EVEN_DATAPATH}/${sep}/*${FILEAPPELLATION})
                        #convert from an array to a... list? ida know. bash stuff.
                        EVEN_FILES=`lst_select.py -C ${cal} --ra=${LST} ${EVEN_FILES[@]}`
                    else
                        EVEN_FILES=${EVEN_DATAPATH}/${sep}/${EVEN_GLOB}.${FILEAPPELLATION}
                    fi

                    #check for odd glob before using lst_select
                    if [[ -z $ODD_GLOB ]]; then
                         ODD_FILES=(${ODD_DATAPATH}/${sep}/*${FILEAPPELLATION})
                         ODD_FILES=`lst_select.py -C ${cal} --ra=${LST} ${ODD_FILES[@]}`
                    else
                         ODD_FILES=${ODD_DATAPATH}/${sep}/${ODD_GLOB}.${FILEAPPELLATION}
                    fi

                    test -e ${sepdir} || mkdir ${sepdir}
                    LOGFILE=`pwd`/${PREFIX}/${chan}_${pol}_${sep}.log
                    echo This is mk_pspec.sh:  |tee -a  ${LOGFILE}
                    echo -e '\t' recording to: ${LOGFILE} | tee -a ${LOGFILE}
                    echo -e '\t' experiment: ${PREFIX}|tee -a ${LOGFILE}
                    echo -e '\t' channels: ${chan}|tee -a ${LOGFILE}
                    echo -e '\t' polarization: ${pol}|tee -a ${LOGFILE}
                    echo -e '\t' separation: ${sep}|tee -a ${LOGFILE}
                    echo -e '\t' date: `date` | tee -a ${LOGFILE}

                    #ANTS=`grid2ant.py -C ${cal} --seps="${sep}"`

                    echo Beginning pspec calculation | tee -a ${LOGFILE}
                    # If plotting covariances
                    if [ $PLOT == True ]; then 
                    echo python ${SCRIPTSDIR}/pspec_cov_v003.py -C ${cal} \
                         -b ${NBOOT} -a ${ANTS} -c ${chan} -p ${pol}\
                          --window=${WINDOW}  ${NOPROJ} --output=${sepdir} --rmbls=${RMBLS} --plot \
                          ${noise_append[${chan}]} \
                           ${EVEN_FILES} ${ODD_FILES} ${OPTIONS}
                    python ${SCRIPTSDIR}/pspec_cov_v003.py -C ${cal} -b ${NBOOT} \
                        -a ${ANTS} -c ${chan} -p ${pol} --window=${WINDOW} \
                          ${NOPROJ} --output=${sepdir} --rmbls=${RMBLS} --plot \
                          ${noise_append[${chan}]} \
                          ${EVEN_FILES} ${ODD_FILES} ${OPTIONS} \
                         | tee -a ${LOGFILE}
                    fi
                    
                    # If not plotting covariances
                    if [ $PLOT == False ]; then
                    echo python ${SCRIPTSDIR}/pspec_cov_v003.py -C ${cal} \
                         -b ${NBOOT} -a ${ANTS} -c ${chan} -p ${pol}\
                          --window=${WINDOW}  ${NOPROJ} --output=${sepdir} --rmbls=${RMBLS} \
                          ${noise_append[${chan}]} \
                           ${EVEN_FILES} ${ODD_FILES} ${OPTIONS}
                    python ${SCRIPTSDIR}/pspec_cov_v003.py -C ${cal} -b ${NBOOT} \
                        -a ${ANTS} -c ${chan} -p ${pol} --window=${WINDOW} \
                          ${NOPROJ} --output=${sepdir} --rmbls=${RMBLS} \
                          ${noise_append[${chan}]} \
                          ${EVEN_FILES} ${ODD_FILES} ${OPTIONS} \
                         | tee -a ${LOGFILE}
                    fi
                if [ $? -ne 0 ] 
                then
                exit
                fi
                fi 
                if [ $BOOT == True ]; then
                echo Beginning bootstrap: `date` | tee -a ${LOGFILE} 
                if [ $PLOT == True ]; then
                ${SCRIPTSDIR}/pspec_cov_boot_v002_pI.py ${sepdir}/pspec_boot*npz  | tee -a ${LOGFILE} 
                else
                ${SCRIPTSDIR}/pspec_cov_boot_v002_pI.py ${sepdir}/pspec_boot*npz  | tee -a ${LOGFILE}
                fi
                fi
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






if [ $KPKPLOT == True ]; then
echo averaging power spectra for pols/channels
for chan in $chans; do
    chandir=${PREFIX}/${chan}
    for pol in $pols; do
        LOGFILE=`pwd`/${PREFIX}/${chan}_${pol}.log
        echo "Generating plots for ${chan}: ${pol}" | tee -a  ${LOGFILE}
        echo "Multiplying pspec by factor ${scales[${chan}]} for Cov" | tee -a  ${LOGFILE}
        poldir=${chandir}/${pol}
        #PLOT
        ${SCRIPTSDIR}/plot_pk_k3pk_mjk.py ${poldir}/*/pspec.npz --cov=${scales[$chan]}  ${nocov} --outpath=${PREFIX}_${chan}_${pol} | tee -a  ${LOGFILE}
        #mv pspec_pk_k3pk.npz pspec.npz
        #mv pspec.png pspec_${PREFIX}_${chan}_${pol}.png
        #mv posterior.png posterior_${PREFIX}_${chan}_${pol}.png
        #mv posterior.txt posterior_${PREFIX}_${chan}_${pol}.txt
        mv pspec_${PREFIX}_${chan}_${pol}.png ${poldir}/
        mv posterior_${PREFIX}_${chan}_${pol}.png ${poldir}/
        mv pspec_${PREFIX}_${chan}_${pol}.npz ${poldir}/
        mv posterior_${PREFIX}_${chan}_${pol}.txt ${poldir}/
    done
done
fi
)
