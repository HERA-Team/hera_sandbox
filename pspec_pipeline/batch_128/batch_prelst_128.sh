#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -j y
#$ -l h_vmem=5G
#$ -N fg_sub
PATH=/home/saulkohn/src/anaconda/bin/:$PATH
source activate PAPER
ARGS=`pull_args.py $*`
CALFILE=psa6622_v003
PATH2CAPO='/home/saulkohn/ReposForCanopy/capo'
echo Processing: $ARGS

for FILE in $ARGS; do
#    echo mk_stokes.py --stokes=I ${FILE}c
#    mk_stokes.py --stokes=I ${FILE}c

    echo ${PATH2CAPO}/pspec_pipeline/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}
    ${PATH2CAPO}/pspec_pipeline/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}

    echo xrfi_simple.py -n 3 ${FILE}B
    xrfi_simple.py -n 3 ${FILE}B && rm ${FILE}B/* && rmdir ${FILE}B
done

echo DONE
