#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -j y
#$ -l h_vmem=800M
#$ -q all.q
#$ -N fg_sub
#$ -t 1:100
PATH=/home/saulkohn/src/anaconda/bin/:$PATH
source activate PAPER
ARGS=`~/scripts/pull_args.py $*`
CALFILE=psa6622_v003
echo Processing: $ARGS

for FILE in $ARGS; do
    echo --------------------
    echo $FILE

#    echo correct_psa6240_omni.py ${FILE}
#    correct_psa6240_omni.py ${FILE}

#    echo mk_stokes.py --stokes=I ${FILE}c
#    mk_stokes.py --stokes=I ${FILE}c

    echo /home/saulkohn/ReposForCanopy/capo/pspec_pipeline/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}
    /home/saulkohn/ReposForCanopy/capo/pspec_pipeline/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}

#    echo xtalk3.py ${FILE}cPB
#    xtalk3.py ${FILE}cPB
   
    echo xrfi_simple.py -n 3 ${FILE}B
    python /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/sak_workspace/xrfi_simple.py -n 3 ${FILE}B && rm ${FILE}B/* && rmdir ${FILE}B
done

echo DONE
