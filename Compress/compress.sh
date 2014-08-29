#! /bin/bash
files2proc=$@
NL="\n"
for f in $files2proc; do
    [[ $f == *npz ]] && continue
    hn=`print_hostname.py`
    infile=${hn}:${f}
    outfile=${infile}E
    
    triplet=`get_neighbor.py ${f}`
    tlen=`echo ${triplet} | wc -w`
    [[ `echo ${triplet} | wc -w` == 3 ]] || continue
    
    Rtriplet=""
    for t in $triplet; do
        [[ -e ${t%R} ]] && Rtriplet="${Rtriplet} ${t%R}" || Rtriplet="${Rtriplet} ${t}" 
    done

    record_launch.py ${outfile} -i ${infile} -d '5-DDR'
    LOG="ddr_filter_coarse.py -a all -p xx,xy,yx,yy --maxbl=301 --clean=1e-3 --nsections=20 $Rtriplet"${NL}
    LOG=${LOG}`date`${NL}
    stdout=$(ddr_filter_coarse.py -a all -p xx,xy,yx,yy --maxbl=300 --clean=1e-3 --nsections=20 ${Rtriplet} 2>&1)
    PID=$!
    STATUS=$?
    LOG=${LOG}${stdout}${NL}
    if [[ $STATUS -eq 0 ]]; then
        for suffix in "D" "E" "F"; do
            test -e ${f}${suffix} && add_file.py ${outfile%E}${suffix} -i ${infile}
        done
        record_completion.py ${outfile} --log="${LOG}"
    else
        rm -r ${f}[DEF]
        record_failure.py ${outfile} --log="${LOG}"
    fi
done

for f in $files2proc; do
    [[ -e ${f}[DEF] ]] && rm -r ${f}
done
