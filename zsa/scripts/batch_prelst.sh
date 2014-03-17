#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=4G
#$ -l paper
ARGS=`pull_args.py $*`
CALFILE=psa6240_v003
echo Processing: $ARGS

for FILE in $ARGS; do
    echo --------------------
    echo $FILE
    #FILE=`python -c "print '${FILE}'[:-1]"`
    ##xrfi_simple.py --dt=3 -c 74,152,153,168 $FILE
    ##xtalk3.py ${FILE}R
#    echo apply_cal.py -C ${CALFILE} ${FILE}
#    apply_cal.py -C ${CALFILE} ${FILE}
    #echo combine_pols.py -p xx,yy ${FILE}C
    #combine_pols.py -p xx,yy ${FILE}C
    ##combine_pols.py -p xy,yx ${FILE}C
#    echo mk_stokes.py --stokes=I ${FILE}C
#    mk_stokes.py --stokes=I ${FILE}C
    #pspec_combine_bls.py -a cross ${FILE} -p xx
#    echo /usr/global/paper/capo/pspec_pipeline/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}CP
#    /usr/global/paper/capo/pspec_pipeline/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}CP
    #/usr/global/paper/capo/pspec_pipeline/pspec_prep.py -C psa898_v003 -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=100 ${FILE}CPS
    #xtalk3.py ${FILE}SB
    #xrfi_simple.py -n 3 -c 0_35,176_202 ${FILE}CPBx
    echo xrfi_simple.py -n 3 ${FILE}CPBx
    xrfi_simple.py -n 3 ${FILE}CPBx
done

echo DONE
