#$ -S /bin/bash
#$ -j y
#$-o grid_output
ARGS=`pull_args.py $*`
CORRECT=correct_pgb600_v005.py

for FILE in $ARGS; do
    #$CORRECT $FILE
    #combine_freqs.py -n 1024 -u ${FILE}c
    #apply_bp.py ${FILE}
    /data1/paper/arp/scripts/xrfi_simple.py --dt=2.5 --df=4 --combine -t 80 ${FILE}
    xtalk3.py ${FILE}R
done
