#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=4G
#$ -l paper
ARGS=`pull_args.py $*`
echo Processing: $ARGS

for FILE in $ARGS; do
    echo --------------------
    mk_stokes.py --stokes=I,Q,U,V ${FILE}
    pspec_combine_bls.py -a cross ${FILE}P
    #/usr/global/paper/capo/pspec_pipeline/pspec_prep.py -C psa898_v003 -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}CPS
    xtalk3.py ${FILE}PS
    xrfi_simple.py -n 4 -c 0_35,176_202 ${FILE}PSx
done

echo DONE
