#! /bin/bash

function branch () 
{
    test -e $1 || mkdir $1
}
function format_jid ()
{
    x=${1##Your job-array}
    x=${x%%.*}
    echo ${x}
}

(
    #Read in the config file
    if [[ $* == "" ]]; then 
        echo "I don't work without a config file."
        exit
    fi
    
    echo "using config ${*}"
    . $*
    #Build directory structure
    echo "Creating directory structure for power spectrum ${PREFIX}"
    branch ${PREFIX}
    for chan in $chans; do
        chandir=${PREFIX}/${chan}
        branch ${chandir}
        for pol in $pols; do
            poldir=${chandir}/${pol}
            branch ${poldir}
            for sep in ${seps}; do
                sepdir=${poldir}/${sep}
                branch ${sepdir}
            done
        done
    done
    
    #fire off jobs.
    jobs=$(ls -d ${PREFIX}/*/*/*)
    t=$(echo ${jobs} | wc -w)
    echo "Submitting ${t} jobs..."
    jid=`qsub -t 1:${t} -q test.q single_job_pspec.sh ${jobs}`
    jid=$(format_jid "${jid}")
    echo "${t} jobs submitted! Waiting for completion..."
    until [[ $(qstat | grep ${jid}) == "" ]]; do
        wait
    done
    
    #average all bootstraps and clean up.
    echo "Jobs done..."
    echo "Averaging power spectra for each pol/channel combination."
    for chan in $chans; do
        chandir=${PREFIX}/${chan}
        for pol in $pols; do
            echo "Generating plots for ${chan}: ${pol}"
            poldir=${chandir}/${pol}
            ${SCRIPTSDIR}/pspec_plot_pk_k3pk.py ${poldir}/*/pspec.npz
            mv pspec_pk_k3pk.npz ${poldir}/*
        done
    done
    echo "DONE!"
)
