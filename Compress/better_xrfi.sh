#! /bin/bash

for f in $@; do
    hn=`hostname`
    infile=${hn}:${f}
    outfile=${infile}R
    
    triplet=`get_neighbor.py ${f}`
    if [ `echo ${triplet} | wc -w` != 3 ]; then
        continue
    fi
    
    if [ ! -f ${f}R ]; then
        #echo record_launch.py ${outfile} -i ${infile} -d '4-XRFI'
        record_launch.py ${outfile} -i ${infile} -d '4-XRFI'
        #ddr_filter_coarse.py -a 1 -p xx,xy,yx,yy --clean=1e-3 --maxbl=300 --output=ddr --invert $triplet
        #xrfi_simple.py -a 1 --combine -t 80 -n 5 ${f}E --to_npz=${f}E.npz
        #xrfi_simple.py -a all --combine -t 80 ${f} --from_npz=${f}E.npz
        dd if=/dev/urandom of=${f}R bs=16 count=1 &> /dev/null
        if [[ $? ]]; then
            #echo add_file.py ${outfile} -i ${infile}
            add_file.py ${outfile} -i ${infile}
            #echo record_completion.py ${outfile}
            record_completion.py ${outfile}
            rm -r ${f}
        else
            echo "DO SOMETHING!"
        fi
    fi
done
