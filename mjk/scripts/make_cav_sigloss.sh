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
               EVEN_FILES=${EVEN_DATAPATH}${sep}/*242.[3456]*uvGL
               ODD_FILES=${ODD_DATAPATH}${sep}/*243.[3456]*uvGL
               test -e ${sepdir} || mkdir ${sepdir}
               LOGFILE=`pwd`/${PREFIX}/${chan}_${pol}_${sep}.log
               echo this is make_sigloss.sh with  |tee  ${LOGFILE}
               echo experiment: ${PREFIX}|tee -a ${LOGFILE}
               echo channels: ${chan}|tee -a ${LOGFILE}
               echo polarization: ${pol}|tee -a ${LOGFILE}
               echo separation: ${sep}|tee -a ${LOGFILE}
               echo `date` | tee -a ${LOGFILE} 
               ANTS='cross'
                   
               echo python ${SCRIPTSDIR}/pspec_cav_v002_sigloss.py \
                    --window=${WINDOW} -p ${pol} -c ${chan} --band=${band}\
                     -b ${NBOOT} -C ${cal} -i ${inject}\
                     -a ${ANTS} --output=${ijectdir} --auto \
                     ${EVEN_FILES} ${ODD_FILES}

               python ${SCRIPTSDIR}/pspec_cav_v002_sigloss.py \
                    --window=${WINDOW} -p ${pol} -c ${chan} --band=${band}\
                     -b ${NBOOT} -C ${cal} -i ${inject} \
                     -a ${ANTS} --output=${injectdir} --auto \
                     ${EVEN_FILES} ${ODD_FILES} \
                     |tee -a ${LOGFILE}
            done
        done
    done
done

echo waiting on `python -c "print len('${PIDS}'.split())"` sigloss threads ${PIDS} 
wait $PIDS


echo "Creating Sigloss Plot"

for chan in $chans; do
    chandir=${PREFIX}/${chan}
    for pol in ${pols}; do
        poldir=${chandir}/${pol}
        python ${SCRIPTSDIR}/plot_sigloss.py ${poldir}/*/inject_*
        mv sigloss.png sigloss_${PREFIX}_${chan}_${pol}.png
        cp sigloss_${PREFIX}_${chan}_${pol}.png ${poldir}/
    done
done    

    )
