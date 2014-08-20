#! /bin/bash
#send files to the sill in a nice orderly pdb friendly way
infile=$1
outfile=$2

if ssh ${still} test -e ${outfile##*:}; then
    continue
else
    record_launch.py ${outfile} -i ${infile} -d '2-RSYNC'
    echo Moving: ${infile}
    LOG="Moving file from pot to still node for processing\n"
    LOG=${LOG}"scp -r -c arcfour256 ${infile} ${outfile}${NL}"
    LOG=${LOG}`date`${NL}
    LOG=${LOG}$(scp -r -c arcfour256 ${infile} ${outfile} 2>&1)
    STATUS=$?
    if [[ $STATUS -eq 0 ]]; then
        ssh ${still} "add_file.py ${outfile} -i ${infile}"
        record_completion.py ${outfile} --log="${LOG}"
    else
        ssh ${still} "rm -r ${outfile##*:}"
        record_failure.py ${outfile} --log="${LOG}"
    fi


