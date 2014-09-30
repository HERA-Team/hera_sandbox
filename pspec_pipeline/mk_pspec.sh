#! /bin/bash
export PYTHONPATH='.':$PYTHONPATH
#PREFIX="OneDayFG"
#
##chans=`python -c "print ' '.join(['%d_%d'%(i,i+39) for i in range(10,150,1)])"`
#pols='I Q U V'
#seps='0_16 1_16 0_17'
#chans='110_149'
#RA="1:01_9:00"
#NBOOT=20
#
##DATAPATH=fringe_hor_v006
#SCRIPTSDIR=~/src/capo/pspec_pipeline
#cal="psa898_v003"
#PWD=`pwd`
#DATAPATH="${PWD}/typical_day/*FRXS"
(
echo using config $*
. $*
#set defaults to parameters which might not be set
if [[ ! -n ${GAIN} ]]; then export GAIN="0.3"; fi
if [[ ! -n ${WINDOW} ]]; then export WINDOW="blackman-harris"; fi
#for posterity Print the cal file we are using

pywhich $cal

threadcount=`python -c "c=map(len,['${pols}'.split(),'${chans}'.split(),'${seps}'.split()]);print c[0]*c[1]*c[2]"`
echo Running $threadcount pspecs
PIDS=""

FILES=`lst_select.py -C ${cal} --ra=${RA} ${DATAPATH}`
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

                test -e ${sepdir} || mkdir ${sepdir}
                LOGFILE=`pwd`/${PREFIX}/${chan}_${pol}_${sep}.log
                echo this is mk_pspec.sh with  |tee  ${LOGFILE}
                echo experiment: ${PREFIX}|tee -a ${LOGFILE}
                echo channels: ${chan}|tee -a ${LOGFILE}
                echo polarization: ${pol}|tee -a ${LOGFILE}
                echo separation: ${sep}|tee -a ${LOGFILE}
                echo gain: ${GAIN} | tee -a ${LOGFILE}
                echo `date` | tee -a ${LOGFILE}

                ANTS=`grid2ant.py -C ${cal} --seps="${sep}"`
                echo ${SCRIPTSDIR}/pspec_redmult_cov_gps.py -C ${cal} -b ${NBOOT} \
                    -a ${ANTS} -c ${chan} -p ${pol} --window=${WINDOW} \
                    --gain=${GAIN} --output=${sepdir} ${NOPROJ} ${FILES} 

                ${SCRIPTSDIR}/pspec_redmult_cov_gps.py -C ${cal} -b ${NBOOT} \
                    -a ${ANTS} -c ${chan} -p ${pol} --window=${WINDOW} \
                    --gain=${GAIN} --output=${sepdir} ${NOPROJ} ${FILES} \
                | tee -a ${LOGFILE}  
                echo beginning bootstrap: `date` | tee -a ${LOGFILE} 
                ${SCRIPTSDIR}/pspec_pk_k3pk_boot.py ${sepdir}/pspec_boot*npz | tee -a ${LOGFILE} 
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
        ${SCRIPTSDIR}/pspec_plot_pk_k3pk.py ${poldir}/*/pspec.npz
        mv pspec_pk_k3pk.npz pspec_${PREFIX}_${chan}_${pol}.npz
        mv pspec.png pspec_${PREFIX}_${chan}_${pol}.png
        cp  pspec_${PREFIX}_${chan}_${pol}.png ${poldir}/
        cp pspec_${PREFIX}_${chan}_${pol}.npz ${poldir}/
    done
done
)
