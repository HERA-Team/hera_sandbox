#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=6G
#$ -l paper
#$ -N fg_sub

ARGS=`pull_args.py $*`
CALFILE=psa6622_v003 #Aaron's quick calfile
echo Processing: $ARGS

for FILE in $ARGS; do
    echo --------------------
    echo $FILE

    if [ -d ${FILE}B ]; then 
        echo ${FILE}B exists... skipping PSPEC_PREP...
    else 
        echo ~/capo/ctc/scripts/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}
        ~/capo/ctc/scripts/pspec_prep.py -C ${CALFILE} -a cross --window=blackman-harris --nogain --nophs --clean=1e-9 --horizon=15 ${FILE}
    fi

    if [ -d ${FILE}BR ]; then
        echo ${FILE}BR exists... skipping XRFI_SIMPLE...
    else
        echo xrfi_simple.py -n 4 ${FILE}B
        xrfi_simple.py -n 4 ${FILE}B #&& rm ${FILE}B/* && rmdir ${FILE}B
    fi

done
echo DONE


