#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=4G
#$ -l paper
ARGS=`pull_args.py $*`
echo batch_lst_pipeline.sh $ARGS

for FILE in $ARGS; do
    echo --------------------
    xrfi_simple.py --dt=3 -c 74,152,153,168 $FILE
    xtalk3.py ${FILE}R
    apply_cal.py -C psa898_v003 ${FILE}Rx
    /usr/global/paper/capo/pspec_pipeline/pspec_prep.py -C psa898_v003 -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}RxC
    xrfi_simple.py -n 3 -c 0_35,176_202 ${FILE}RxCB
    xtalk3.py ${FILE}RxCBR
    combine_pols.py -p xx,yy ${FILE}RxCBRx
    pspec_combine_bls.py -a cross ${FILE}RxCBRxP
done

echo DONE
