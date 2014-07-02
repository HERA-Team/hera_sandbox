#! /bin/bash

SCRATCH=/data
#moved to correct_and_xrfi.sh
RFI_CHANS="0_65,377_388,510,770,840,852,913,921_922,932_934,942_1023"
Nproc=8

function manytimes {
    #hacky rewrite of qsub ---      
    n=0
    PIDS=""
    while [[ $n -lt $Nproc ]]; do
        export n
        export Nproc
        _f=`file2core.py $files`
        $@ $_f &
        PIDS="${PIDS} "$! 
        n=$((n+1))
    done
    wait $PIDS
}

#3-XRFI
#rm
(
    files=`ls -d ${SCRATCH}/zen*uv`
    export files
    manytimes /home/obs/Compress/correct_and_xrfi.sh 
    wait $! 
)
#4-XRFI
#rm
(
    files=`ls -d ${SCRATCH}/zen*cR`
    export files
    manytimes /home/obs/Compress/better_xrfi.sh 
    wait $! 
)
#5-DDR
#rm
(
    files=`ls -d ${SCRATCH}/zen*cR*`
    export files
    manytimes /home/obs/Compress/compress.sh 
    wait $! 
)

#6-RSYNC 
/home/obs/Compress/clean_up.sh
