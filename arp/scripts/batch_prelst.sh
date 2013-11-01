#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=4G
#$ -l paper
ARGS=`pull_args.py $*`
echo Processing: $ARGS

for FILE in $ARGS; do
    echo --------------------
    echo $FILE
    FILE=`python -c "print '${FILE}'[:-1]"`
    ##xrfi_simple.py --dt=3 -c 74,152,153,168 $FILE
    ##xtalk3.py ${FILE}R
    #apply_cal.py -C psa898_v003 ${FILE}
    ##combine_pols.py -p xx,yy ${FILE}C
    ##combine_pols.py -p xy,yx ${FILE}C
    mk_stokes.py --stokes=I ${FILE}C
    pspec_combine_bls.py -a cross ${FILE}CP
    /usr/global/paper/capo/pspec_pipeline/pspec_prep.py -C psa898_v003 -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}CPS
    #/usr/global/paper/capo/pspec_pipeline/pspec_prep.py -C psa898_v003 -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=100 ${FILE}CPS
    xtalk3.py ${FILE}CPSB
    xrfi_simple.py -n 3 -c 0_35,176_202 ${FILE}CPSBx
done

echo DONE
