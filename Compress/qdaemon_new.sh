#! /bin/bash

function send_file_to_stillhost() {
    infile=$1
    outfile=$2
    send_file_to_stillhost.sh $infile $outfile &
    PID=$!
    wait_for_still_process_start.py "2-RSYNC" $infile
    record_pid.py --pid=${PID} $infile
}

function remote_launch() {
    COMMAND=$1
    HOST=$2
    FILE=$3
    TASK=$4
    PID=(ssh -f $HOST "nohup $COMMAND $FILE &; echo $!")
    wait_for_still_process_start.py $TASK $FILE #TODO check that this is the right file to be waiting for. output vs basefile
    record_pid.py --pid=${PID} $FILE
}  
#newline character for legible formatting
#NL=$'\n'
#where is scratch space on still machines?
SCRATCH=/data
MaxFiles=288 #the maximum number of files I ever want to send to the still at once.

#Read in the list of available compute nodes, store them to a variable.
declare -a ComputeNodes
readarray -t ComputeNodes < /home/obs/groups/compute_nodes
Nnodes=${#ComputeNodes[@]}

while :
do
    #first check that there are non edge files to compress
    UncompressedFiles=(`get_files_to_distill.py ${MaxFiles}`)
    if [[ ${#UncompressedFiles[*]} -eq 0 ]]; then continue fi  #this step prevents us from looping forever on edge files
    #what files do I need to compress? put them in a bash array. include the edges
    Observations2Compress=(`get_files_to_distill.py --include_edges ${MaxFiles}`)
    Nfiles=${#Observations2Compress[@]}
     
    #only work if I need to.
    if [[ $Nfiles -gt 0 ]]; then
        echo "Found ${Nfiles}, beginning compression sequence"
        DONT_STOP='true'
        #until I get the stop condition...
        while $DONT_STOP;  do
            DONT_STOP=`check_still_stop_condition.py ${Observations2Compress[@]}`
            AVAILABLE_SLOTS=`get_available_slots.py` 
            #storage=files table, if a basefile is on a host, counts as one slot (no matter how many)
            #Observations2Compress=`get_observations_to_distill.py $AVAILABLE_SLOTS`
            #Observations2Compress=${Observations2Compress}`get_current_distilling.py`
            #assign files by first random + sticky.
            #loop through all the _observations_  
            for ((i=0;i<$Nfiles;i++)); do  #loop instead Observations2Compress
                #which still does this live on?
                #recalculated for 2-RSYNC  
                #for each observation, loop through all the most recent files that have been generated and do any work 
                # required
                Files2Compress=`get_latest_completed_files.py $Observations2Compress[$i]}`
                for File2Compress in $Files2Compress; do #TODO r/${Files2Compress[$i]}/$File2Compress/g
                    still=`get_host.py ${Files2Compress[$i]}`
                    #still=${File2Compress%%:*}
                    #decide what the next step is
                    ActionItem=`get_next_step.py ${Files2Compress[$i]}`
                    #pass the proper command
                    case $ActionItem in
                        '2-RSYNC') #pot?-->still?
                            #assign a still
                            places2go=`file2still.py ${i} ${Nfiles} ${Nnodes}` #returns a list of still numbers
                            for si in $places2go; do
                                still=${ComputeNodes[$si]}
                                infile=${Files2Compress[$i]}
                                outfile=${still}:${SCRATCH}/${infile##*/}
                                #a suggested way to set up the file distribution. --DCJ
                                #send Files2Compress[$i] and its neighbors
                                send_file_to_stillnode ${infile} ${outfile} #note this is a bash wrapper around a shell script of the same name
                                Neighbors=`get_neighbors ${infile}`
                                for infile in $Neighbors; do
                                    outfile=${si}:${SCRATCH}/${infile##*/}
                                    send_file_to_stillnode $infile $outfile
                                done
                            done
                            ;;
                        '3-XRFI')
                            #fire off first-round xrfi
                            ssh -f ${still} "correct_and_xrfi.sh ${Files2Compress[$i]}"
                            #remote_launch "correct_and_xrfi.sh" ${still} 
                            ;;
                        '4-XRFI')
                            #fire off second-round xrfi
                            ssh -f ${still} "better_xrfi.sh ${Files2Compress[$i]}"
                            ;;
                        '5-DDR')
                            #fire off compress job
                            ssh -f ${still} "compress.sh ${Files2Compress[$i]}"
                            ;;
                        '6-RSYNC')
                            
                            ;;
                        '7-RSYNC')
                            #Send to USA|ASU
                            ;;
                        'RESET')
                            ;;
                        'ERROR')
                            ;;
                        'NULL')
                            continue
                            ;;
                    esac
                done
            done
        done
        #once the stop condition is satisfied, I can now clean up.
        for still in "${ComputeNodes[@]}"; do
            ssh -f ${still} clean_up.sh
        done
    else
        #if I'm waiting on files to arrive, pass this handy message.
        now=`date -u`
        echo "Waiting for new data (${now})"
        sleep 10
    fi
done
