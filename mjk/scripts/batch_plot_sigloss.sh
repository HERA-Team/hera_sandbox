#! /bin/bash

export pythonpath='.':$pythonpath
(
echo using config $*
. $*
#set defaults to parameters which might not be set
if [[ ! -n ${window} ]]; then export window="none"; fi
#for posterity print the cal file we are using

pywhich $cal
echo "Creating Sigloss Plot"

for chan in $chans; do
    chandir=${PREFIX}/${chan}
    for pol in ${pols}; do
        poldir=${chandir}/${pol}
        LOGFILE=`pwd`/${PREFIX}/${chan}_${pol}.log
        echo Now plotting sigloss for ${chan}  | tee ${LOGFILE}
        echo Polarization: ${pol} | tee -a ${LOGFILE}
        pspecdir=${PSPEC}/${chan}/${pol}
        echo   python ${SCRIPTSDIR}/plot_sigloss.py ${poldir}/*/inject_* \
            --path=${pspecdir} --pspec=${PSPEC} --pol=${pol} \
             | tee -a ${LOGFILE}
        python ${SCRIPTSDIR}/plot_sigloss.py ${poldir}/*/inject_* \
            --path=${pspecdir} --pspec=${PSPEC} --pol=${pol}\
            | tee -a ${LOGFILE}
        mv sigloss.png sigloss_${PREFIX}_${chan}_${pol}.png
        mv sigloss_${PREFIX}_${chan}_${pol}.png ${poldir}/
    done
done    

)
