#! /bin/bash

SCRATCH=/data
hn=`print_hostname.py`

for f in `ls -d $SCRATCH/*`; do
    if [[ "${f}" == *[DEFz] ]]; then
        infile=${hn}:${f}
        pot=`get_base_host.py ${infile}`
        FinalDir=`get_base_dir.py ${infile}`
        outfile=${FinalDir}/${f##*/}
        if ! ssh ${pot} test -e ${outfile##*:}; then
            record_launch.py ${outfile} -i ${infile} -d '6-RSYNC'
            LOG="scp ${infile} ${outfile}\n"
            LOG=${LOG}`date`"\n"
            LOG=${LOG}$(scp -r -c arcfour256 ${infile} ${outfile} 2>&1)"\n"
            PID=$!
            STATUS=$?
            if [[ $STATUS -eq 0 ]]; then
                ssh ${pot} "add_file.py ${outfile} -i ${infile}"
                record_completion.py ${outfile} --log="${LOG}"
            else
                record_failure.py ${outfile} --log="${LOG}"
            fi
        fi
    fi
   
    delete_order.py ${hn}:${f} '6-RSYNC'
    #^^^ I don't want to do this yet.
    echo "rm -r ${f}"
    rm -r ${f}
done
