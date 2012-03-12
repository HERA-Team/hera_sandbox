#$ -S /bin/bash
#$ -N pgb050_c
#$ -j y
ARGS=`pull_args.py $*`
if ! ls /scratch/paper &> /dev/null; then
    mkdir /scratch/paper
fi
LOCAL=/scratch/paper/pgb050/
if ! ls $LOCAL &> /dev/null; then
    echo Creating local data dir...
    mkdir $LOCAL
#    mv /tmp/lst*uv $LOCAL
fi
for F in $ARGS; do
    echo Clearing old data...
    rm -r ${LOCAL}${F}c
    if ! ls ${LOCAL}${F} &> /dev/null;then
       echo preloading ${F}...  
       scp -rC ${F} ${LOCAL}
    fi
    TARGS=${LOCAL}`python -c "print '${F}'.split('/')[-1]"` ${TARGS}
done
SRC=`python -c "print '${ARGS[0]}'.split('/')[:-1]"`

echo correct_psa113.py $TARGS
correct_psa113.py  $TARGS
echo Completed correction, leaving data on node${SGE_TASK_ID}
