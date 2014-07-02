#! /bin/bash

files2proc=$@

for f in $files2proc; do
    hn=`print_hostname.py`
    infile=${hn}:${f}
    outfile=${infile}R
    
    triplet=`get_neighbor.py ${f}`
    if [ `echo ${triplet} | wc -w` != 3 ]; then
        continue
    fi

    if [ ! -f ${f}R ]; then
        record_launch.py ${outfile} -i ${infile} -d '4-XRFI'
        echo ddr_filter_coarse.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $triplet
        stdout1=$(ddr_filter_coarse.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert ${triplet} 2>&1)
        echo xrfi_simple.py -a 1 --combine -t 80 -n 5 ${f}E --to_npz=${f}E.npz
        stdout2=$(xrfi_simple.py -a 1 --combine -t 80 -n 5 ${f}E --to_npz=${f}E.npz 2>&1)
        if [[ $? -eq 0 ]]; then
            rfifile=${infile}E.npz
            add_file.py ${rfifile} -i ${infile}
        fi
        LOG="xrfi_simple.py -a all --combine -t 80 ${f} --from_npz=${f}E.npz\n"
        LOG=${LOG}`date`"\n"
        LOG=${LOG}$(xrfi_simple.py -a all --combine -t 80 ${f} --from_npz=${f}E.npz 2>&1)"\n"
        STATUS=$?
        PID=$!
        #dd if=/dev/urandom of=${f}R bs=16 count=1 &> /dev/null
        if [[ $STATUS -eq 0 ]]; then
            add_file.py ${outfile} -i ${infile}
            record_completion.py ${outfile} --log="${LOG}"
            #rm -r ${f}
        else
            rm -r ${f}R
            record_failure.py ${outfile} --log="${LOG}"
        fi
    fi
done
#garbage collection
for f in $files2proc; do
    if [ -e ${f}E ]; then
        rm -r ${f}E
    fi
    if [ -e ${f}R ]; then
        rm -r ${f}
    fi
done
