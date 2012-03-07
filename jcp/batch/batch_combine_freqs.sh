#$ -S /bin/bash
#$ -j y
#$-o grid_output
ARGS=`pull_args.py $*`
#echo $ARGS

for FILE in $ARGS; do
    #echo $FILE 
    #$CORRECT $FILE
    combine_freqs.py -n 256 ${FILE}
    #apply_bp.py ${FILE}cm
    #/data1/paper/arp/scripts/xrfi_simple.py --dt=2.5 --df=4 --combine -t 5 ${FILE}cmb
    #xtalk3.py ${FILE}cmbR
done
