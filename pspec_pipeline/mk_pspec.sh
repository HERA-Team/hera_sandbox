#! /bin/bash
PREFIX="TestPlumbing"

#chans=`python -c "print ' '.join(['%d_%d'%(i,i+39) for i in range(10,150,1)])"`
pols='I Q U V'
seps='0_16 1_16 0_17'
chans='110_149'
RA="1:01_9:00"
NBOOT=20

#DATAPATH=fringe_hor_v006
SCRIPTSDIR=~/scripts/
cal="psa898_v003"
PWD=`pwd`
DATAPATH="${PWD}/lst*uv"
PIDS=""

FILES=`${SCRIPTSDIR}/lst_select.py -C ${cal} --ra=${RA} ${DATAPATH}`
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
                echo "   Working on ${pol}: ${sep}"
                cd ${sepdir}
                ${SCRIPTSDIR}/pspec_redmult_cov_gps.py -b ${NBOOT} -a ${sep} -c ${chan} -p ${pol} ${FILES} \
                    && ${SCRIPTSDIR}/pspec_pk_k3pk_boot.py pspec_boot*npz &
                PIDS="${PIDS} "$!
                cd -    
            done
        fi
    done
done

echo waiting on power spectra ${PIDS}
wait $PIDS

echo averaging power spectra for channels
for chan in $chans; do
    chandir=${PREFIX}/${chan}
    for pol in $pols; do
        echo "Generating plots for ${chan}: ${pol}"
        poldir=${chandir}/${pol}
        #PLOT
        ${SCRIPTSDIR}/pspec_plot_pk_k3pk.py ${PREFIX}/${chan}/pspec_*/pspec.npz
        mv pspec_pk_k3pk.npz ${poldir}/
        mv pspec.png pspec_${PREFIX}_${chan}_${pol}.png
        cp  pspec_${PREFIX}_${chan}_${pol}.png ${poldir}/
    done
done
