#! /bin/bash

"""
daemon script overseeing compression.
 1) read trigger file.
 2) make new output directory for compression 
 3) allocate files to each compute host. scp them.
 4) begin compression on each host.

To Do:
  --- remove pol-specfic allocation
  --- incorporate database interface.
  --- Update plumbing step at the end. Should the correlator put files in their final resting place?

DFM 
"""

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
            echo "Wipint scratch from ${HOST}"
            ssh -q -o ConnectTimeout=3 ${HOST} "rm -r ${SCRATCH}/*"
            while read line; do
                if ssh ${HOST} test -e ${SCRATCH}/${line##*} < /dev/null; then
                    echo "Target File(${HOST}:${SCRATCH}/${line##*/}) Exists!"
                else
                    if [[ $line == /data0/* ]]; then
                        scp -rp -c arcfour256 pot0:${line} ${HOST}:${SCRATCH}
                    else
                        scp -rp -c arcfour256 pot1:${line} ${HOST}:${SCRATCH}
                    fi
                fi
                ssh ${HOST} test -e ${SCRATCH}/${line##*/} < /dev/null || echo $line >> ${OUtputDdir}/MIssingFiles.txt
            done < $TRIGGER
            jid1=`qsub -t 1:10 -q ${QUEUE} -o ${OutputDir} -j y correct_and_xrfi.sh`    
            jid1=`format qsub "${jid1}"`
            jid2=`qsub -t 1:10 -q ${QUEUE} -o ${OutputDir} -j y -hold_jid ${jid2} better_xrfi.sh`    
            jid2=`format qsub "${jid12}"`
            jid3=`qsub -t 1:10 -q ${QUEUE} -o ${OutputDir} -j y -hold_jid ${jid3} compress.sh`    
            CurrentJobs[$i]=$jid3
        done
        ssh -q -o ConnectTimeout=3 pot0 "/home/obs/daily_move_pot0.sh ${thisJD}"
        tend=`date -u`
        echo "Tx begin: ${tstart}\nTx end: ${tend}" > ${OutputDir}/TransferLog.txt
    fi
done
