#! /bin/bash

SCRATCH=/data

function manytimes {
    #hacky rewrite of qsub ---      
    n=0
    Nproc=8
    while [[ $n -lt $Nproc ]]; do
        (
            export n
            export Nproc
            _f=`qsub_hack.py $files`
            $@ $_f
        )
        n=$((n+1))
    done
}



files=`ls -d ${SCRATCH}/zen*`
(
    export files
    manytimes echo 
)
#3-XRFI
#4-rm
#5-XRFI
#6-rm
#7-DDR
#8-RSYNC
