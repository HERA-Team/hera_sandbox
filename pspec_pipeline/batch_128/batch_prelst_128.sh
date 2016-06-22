#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -j y
#$ -l h_vmem=800M
#$ -q all.q
#$ -N fg_sub
PATH=/home/jacobsda/src/anaconda/bin/:$PATH
source activate PAPER
ARGS=`~/scripts/pull_args.py $*`
CALFILE=psa6622_v001
echo Processing: $ARGS

for FILE in $ARGS; do
    echo --------------------
    echo $FILE

#    echo correct_psa6240_omni.py ${FILE}
#    correct_psa6240_omni.py ${FILE}

#    echo mk_stokes.py --stokes=I ${FILE}c
#    mk_stokes.py --stokes=I ${FILE}c

    echo /usr/global/paper/capo/zsa/scripts/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}
    ~/src/capo/pspec_pipeline/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}

#    echo xtalk3.py ${FILE}cPB
#    xtalk3.py ${FILE}cPB
   
    echo xrfi_simple.py -n 3 ${FILE}B
    xrfi_simple.py -n 3 ${FILE}B && rm ${FILE}B/* && rmdir ${FILE}B
done

echo DONE
