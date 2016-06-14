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
        test -e $poldir || mkdir $poldir
        for sep in $seps; do 
           sepdir=${poldir}/${sep}
           test -e $sepdir || mkdir $sepdir
           for inject in `python -c "import numpy; print ' '.join(map(str,numpy.logspace(-2,2,40)))"`; do
               injectdir=${sepdir}/inject_${inject}
               test -e ${injectdir} || mkdir ${injectdir}
               echo ${inject}
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
               echo this is make_sigloss.sh with  |tee  ${LOGFILE}
               echo experiment: ${PREFIX}|tee -a ${LOGFILE}
               echo channels: ${chan}|tee -a ${LOGFILE}
               echo polarization: ${pol}|tee -a ${LOGFILE}
               echo separation: ${sep}|tee -a ${LOGFILE}
               echo `date` | tee -a ${LOGFILE} 
                   
               echo python ${SCRIPTSDIR}/pspec_cov_v002_sigloss.py \
                    --window=${WINDOW} -p ${pol} -c ${chan} -b ${NBOOT} \
                     -C ${cal} -i ${inject} -a ${ANTS}\
                     --sep=${sep} --output=${ijectdir} \
                     --frpad=${FRPAD} \
                     ${EVEN_FILES} ${ODD_FILES} 

               python ${SCRIPTSDIR}/pspec_cov_v002_sigloss.py \
                    --window=${WINDOW} -p ${pol} -c ${chan} -b ${NBOOT} \
                     -C ${cal} -i ${inject} -a ${ANTS}\
                     --sep=${sep} --output=${injectdir} \
                     --frpad=${FRPAD} \
                     ${EVEN_FILES} ${ODD_FILES} \
                     |tee -a ${LOGFILE}
            done
        done
    done
done

echo waiting on `python -c "print len('${PIDS}'.split())"` power spectra threads ${PIDS} 
wait $PIDS


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
