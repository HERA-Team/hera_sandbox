#$ -S /bin/bash
#$ -N psa332_b
#$ -j y
#$ -o grid_output/
#$ -V 
#$ -cwd
ARGS=`pull_args.py $*`
#if ! ls /scratch/paper &> /dev/null; then
#    mkdir /scratch/paper
#fi
#LOCAL=/scratch/paper/pgb113/
#if ! ls $LOCAL &> /dev/null; then
#    echo Creating local data dir...
#    mkdir $LOCAL
##    mv /tmp/lst*uv $LOCAL
#fi
#for F in $ARGS; do
##    echo Clearing old data...
##    rm -r ${LOCAL}${F}
#    #if ! ls ${LOCAL}${F} &> /dev/null;then
#    echo preloading ${F}...  
#    scp -r ${F} ${LOCAL}
#    #fi
#    TARGS=${LOCAL}`python -c "print '${F}'.split('/')[-1]"` ${TARGS}
#done
#SRC=`python -c "print '${ARGS[0]}'.split('/')[:-1]"`

echo apply_bp.py -l null -s 6.16327e-8 $ARGS
apply_bp.py  -l null -s 6.16327e-8 $ARGS
#echo Completed bandpass application, leaving data on node${SGE_TASK_ID}
