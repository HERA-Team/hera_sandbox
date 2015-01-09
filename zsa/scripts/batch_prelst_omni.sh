#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=5G
#$ -l paper
ARGS=`pull_args.py $*`
CALFILE=psa6240_v003
echo Processing: $ARGS

for FILE in $ARGS; do
    echo --------------------
    echo $FILE

    echo correct_psa6240_omni.py ${FILE}
    correct_psa6240_omni.py ${FILE}

    echo mk_stokes.py --stokes=I ${FILE}c
    mk_stokes.py --stokes=I ${FILE}c

    echo /usr/global/paper/capo/zsa/scripts/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}cP
    /usr/global/paper/capo/zsa/scripts/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}cP

    echo xtalk3.py ${FILE}cPB
    xtalk3.py ${FILE}cPB

    echo xrfi_simple.py -n 3 ${FILE}cPBx
    xrfi_simple.py -n 3 ${FILE}cPBx
done

echo DONE
