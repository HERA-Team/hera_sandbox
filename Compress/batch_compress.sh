#! /bin/bash

SCRATCH=/data
RFI_CHANS="0_65,377_388,510,770,840,852,913,921_922,932_934,942_1023"
Nproc=8

function manytimes {
    #hacky rewrite of qsub ---      
    n=0
    PIDS=""
    while [[ $n -lt $Nproc ]]; do
        export n
        export Nproc
        _f=`qsub_hack.py $files`
        $@ $_f &
        PIDS="${PIDS} "$! 
        n=$((n+1))
    done
    wait $PIDS
}

#3-XRFI
#rm
(
    files=`ls -d ${SCRATCH}/zen*`
    export files
    manytimes /home/obs/Compress/correct_and_xrfi.sh 
)
#4-XRFI
#rm
(
    files=`ls -d ${SCRATCH}/zen*`
    export files
    manytimes /home/obs/Compress/better_xrfi.sh 
)
#5-DDR
#rm
(
    files=`ls -d ${SCRATCH}/zen*`
    export files
    manytimes /home/obs/Compress/compress.sh 
)

#6-RSYNC 
files=`ls -d ${SCRATCH}/zen*E*`
hn=`hostname`
for f in $files; do
    infile=${hn}:${f}
    pot=`get_base_host.py ${infile}`
    FinalDir=`get_base_dir.py ${infile}`
    outfile=${FinalDir}/${f##*/}
    record_launch.py ${outfile} -i ${infile} -d '6-RSYNC'
    echo scp ${infile} ${outfile}
    scp ${infile} ${outfile}
    if [[ $? ]]; then
        ssh ${pot} "add_file.py ${outfile} -i ${infile}"
        ssh ${pot} "record_completion.py ${outfile}"
    else
        rm -r ${f}
    fi
done
#if anything is left in scratch, kill it.
rm -r ${SCRATCH}/zen*









