#! /bin/bash

#daemon script overseeing compression.
# 1) read in the orders that the correlator sends.
#   w output directory for compression 
# 3) allocate files to each compute host. scp them.
# 4) begin compression on each host.
#
#To Do:
#  --- Update plumbing step at the end. Should the correlator put files in their final resting place?
#
#DFM 

#TriggerRegexp='/home/obs/Compress/Trigger*' #TEST
ARXIVDIR=/home/obs/logs #TEST
SCRATCH=/data
NL="\n"

declare -a ComputeNodes
readarray -t ComputeNodes < /home/obs/groups/compute_nodes
Nnodes=${#ComputeNodes[@]}

while : 
do
    #if the still is running, hold.
    #only generate list of recent orders once the still is finished running.
    while ! `is_still_idle.py`; do
        sleep 10
    done
    
    #read in list of files to compress ---- store them in an array. 
    Trigger=`get_recent_orders.py`
    
    #declare -a files2compress
    files2compress=(${Trigger// / })
    #files2compress=(`get_files_to_distill.py 300`)
    Nfiles=${#files2compress[*]}

    #only try to run if there are files to compress.
    if [[ $Nfiles -gt 0 ]]; then
        #check pdb for incomplete histories.
        echo "${Nfiles} uncompressed files detected, beginning transfer"

        #generate a directory for log messages
        ustart=`date -u "+%F_%H.%M.%S"`
        OutputDir=${ARXIVDIR}/${ustart}
        test -e ${OutputDir} || mkdir ${OutputDir}
        
        for still in ${ComputeNodes[@]}; do
            for ((i=0;i<$Nfiles;i++)); do
                allstills=`file2still.py ${i} ${Nnodes} ${Nfiles}`
                infile=${files2compress[$i]}
                outfile=${still}:${SCRATCH}/${infile##*/}
                for si in $allstills; do
                    if [[ ${ComputeNodes[$si]} == $still ]]; then
                        #2-RSYNC
                        if ssh ${still} test -e ${outfile##*:}; then
                            continue
                        else
                            record_launch.py ${outfile} -i ${infile} -d '2-RSYNC' 
                            echo Moving: ${infile}
                            LOG="Moving file from pot to still node for processing\n"
                            LOG=${LOG}"scp -r -c arcfour256 ${infile} ${outfile}${NL}"
                            LOG=${LOG}`date`${NL}
                            LOG=${LOG}$(scp -r -c arcfour256 ${infile} ${outfile} 2>&1)
                            PID=$!
                            STATUS=$?
                            if [[ $STATUS -eq 0 ]]; then
                                ssh ${still} "add_file.py ${outfile} -i ${infile}" 
                                record_completion.py ${outfile} --log="${LOG}"
                            else
                                ssh ${still} "rm -r ${outfile##*:}"
                                record_failure.py ${outfile} --log="${LOG}"
                            fi
                        fi
                    fi
                done
            done
            ssh -f ${still} batch_compress.sh
        done
    else
        now=`date -u`
        echo "Waiting for new data (${now})"
        sleep 10
    fi
done
