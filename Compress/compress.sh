#! /bin/bash

for f in $@; do
    hn=`hostname`
    infile=${hn}:${f}
    outfile=${infile}E
    
    triplet=`get_neighbor.py ${f}`
    Rtriplet=""
    for t in $triplet; do
        [[ -e ${t%R} ]] && Rtriplet="${Rtriplet} ${t%R}" || Rtriplet="${Rtriplet} ${t}" 
    done

    #echo record_launch.py ${outfile} -i ${infile} -d '5-DDR'
    record_launch.py ${outfile} -i ${infile} -d '5-DDR'
    #ddr_filter_coarse.py -a all -p xx,xy,yx,yy --maxbl=300 --clean=1e-3 --nsections=20 $Rtriplet
    dd if=/dev/urandom of=${f}E bs=16 count=1 &> /dev/null
    if [[ $? ]]; then
        #echo add_file.py ${outfile} -i ${infile}
        add_file.py ${outfile} -i ${infile}
        #echo record_completion.py ${outfile}
        record_completion.py ${outfile}
        rm -r ${f}
    else
        echo "DO SOMETHING!"
    fi
done
