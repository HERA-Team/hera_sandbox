#! /bin/bash

#daemon script overseeing compression.
# 1) read in the orders that the correlator sends.
# 2) make new output directory for compression 
# 3) allocate files to each compute host. scp them.
# 4) begin compression on each host.
#
#To Do:
#  --- Update plumbing step at the end. Should the correlator put files in their final resting place?
#
#DFM 

TriggerRegexp='/home/obs/Compress/Trigger*' #TEST
ARXIVDIR=/home/obs/logs #TEST
SCRATCH=/data

declare -a ComputeNodes
readarray -t ComputeNodes < /home/obs/groups/compute_nodes
Nnodes=${#ComputeNodes[@]}

while : 
do
    #read in list of files to compress ---- store them in an array. 
    Trigger=`get_recent_orders.py '1-RSYNC'`
    
    declare -a files2compress
    files2compress=(${Trigger// / })
    Nfiles=${#files2compress[@]}
   
    #only try to run if there are files to compress.
    if [[ $Nfiles -gt 0 ]]; then
        echo "${Nfiles} uncompressed files detected, beginning transfer"
        sleep 1

        #generate a directory for log messages
        ustart=`date -u "+%F_%H.%M.%S"`
        OutputDir=${ARXIVDIR}/${ustart}
        test -e ${OutputDir} || mkdir ${OutputDir}
        
        for still in ${ComputeNodes[@]}; do
            for ((i=0;i<$Nfiles;i++)); do
                allstills=`file2still.py ${i} ${Nnodes} ${Nfiles}`
                infile=${files2compress[$i]}
                outfile=${still}:${SCRATCH}/${infile##*/}
                #pot=`get_host.py ${infile}`
                for si in $allstills; do
                    if [[ ${ComputeNodes[$si]} == $still ]]; then
                        #2-RSYNC
                        ssh ${still} "record_launch.py ${outfile} -i ${infile} -d '2-RSYNC'"
                        scp ${infile} ${outfile}
                        if [[ $? ]]; then
                            ssh ${still} "add_file.py ${outfile} -i ${infile}" 
                            ssh ${still} "record_completion.py ${outfile}" 
                        else
                            echo DO SOMETHING!
                        fi
                    fi
                done
            done
            ssh ${still} batch_compress.sh
        done
    fi
done
