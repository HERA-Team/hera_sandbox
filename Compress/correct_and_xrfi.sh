#! /bin/bash
NL="\n"
RFI_CHANS="0_65,377_388,510,770,840,852,913,921_922,932_934,942_1023"
for f in $@; do
    hn=`print_hostname.py`
    infile=${hn}:${f}
    outfile=${infile}cR
    
    #echo record_launch.py ${outfile} -i ${infile} -d '3-XRFI'
    record_launch.py ${outfile} -i ${infile} -d '3-XRFI'
    LOG="correct_and_XRFI.py -a 1 -t 80 --df=6 -c ${RFI_CHANS} ${f}  2>&1"
    LOG=${LOG}`date`${NL}
    LOG=${LOG}$(correct_and_XRFI.py -a 1 -t 80 --df=6 -c ${RFI_CHANS} ${f} 2>&1)${NL}
    if [[ $? -eq 0 ]]; then 
        #echo add_file.py ${outfile} -i ${infile}
        add_file.py ${outfile} -i ${infile}
        #echo record_completion.py ${outfile}
        record_completion.py ${outfile} --log="${LOG}"
        rm -r ${f}
    else
        rm -r ${f}cR
        record_failure.py ${outfile} --log="${LOG}"
    fi
done
