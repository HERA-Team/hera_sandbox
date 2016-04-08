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
        pspecdir=${PSPEC}/${chan}/${pol}
        python ${SCRIPTSDIR}/plot_sigloss.py ${poldir}/*/inject_* \
            --path=${pspecdir} --pspec=${PSPEC} --pol=${pol}
        mv sigloss.png sigloss_${PREFIX}_${chan}_${pol}.png
        cp sigloss_${PREFIX}_${chan}_${pol}.png ${poldir}/
    done
done    

)
