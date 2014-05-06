#1 /bin/bash

for f in $@; do
    hn=`hostname`
    infile=${hn}:${f}
    outfile=${infile}cR
    
    #echo record_launch.py ${outfile} -i ${infile} -d '3-XRFI'
    record_launch.py ${outfile} -i ${infile} -d '3-XRFI'
    #echo correct_and_XRFI.py -a 1 -t 80 --df=6 -c ${RFI_CHANS} ${f}
    dd if=/dev/urandom of=${f}cR bs=16 count=1 &> /dev/null
    if [[ $? ]]; then 
        #echo add_file.py ${outfile} -i ${infile}
        add_file.py ${outfile} -i ${infile}
        #echo record_completion.py ${outfile}
        record_completion.py ${outfile}
        rm -r ${f}
    else
        rm -r ${f}cR
        record_failure.py ${outfile}
    fi
done
