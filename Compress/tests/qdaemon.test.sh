#! /bin/bash

TriggerRegexp='Trigger*' #TEST
ARXIVDIR=tests/grid_output #TEST
SCRATCH=/scratch

function format_qsub () {
    echo `echo $1 | awk '{print $3}'`
}

POLS=('' 'xx' 'xy' 'yx' 'yy') #TEST
CurrentJobs=('' '' '' '' '') #TEST

while : 
do
    Trigger=`ls -1dr ${TriggerRegexp} 2> /dev/null | head -n 1`
    if [[ $Trigger != "" ]]; then
        thisJD=${Trigger##*_}; thisJD=${thisJD%%.txt}
        echo "trigger ${Trigger} detected, beginning transfer"
        sleep 1
        tstart=`date -u`
        ustart=`date -u "+%F_%H.%M.%S"`
        OutputDir=${ARXIVDIR}/${ustart}
        test -e ${OutputDir} || mkdir ${OutputDir}

        TRIGGER=${OutputDir}/${Trigger##*/}
        echo; echo "Moving trigger $Trigger to archive $TRIGGER"; echo
        for i in {1..4}; do #TEST --- change this to be more general
            HOST=still${i}
            QUEUE=S${i}.q
            until [[ `qstat | grep ${CurrentJobs[$i]}` == "" ]]; do
                sleep 3
            done
        done
        echo "Wipint scratch from ${HOST}"
        ssh -q -o ConnectTimeout=3 ${HOST} "rm -r ${SCRATCH}/*"
        while read line; do
            if ssh ${HOST} test -e ${SCRATCH}/${line##*} < /dev/null; then
                echo "Target File(${HOST}:${SCRATCH}/${line##*/}) Exists!"
            else
                if [[ $line == /data0/* ]]; then
                    echo data0
                else
                    echo !data0
                fi
            fi
        done < $TRIGGER
    fi
done
