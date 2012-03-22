#$ -S /bin/bash
#$ -j y
#$ -o /data1/paper/jacobsda/stokes/pgb050/grid_output/
#$ -N pgb050_pipe
ARGS=`pull_args.py -w $*`
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
    if ! ls ${LOCAL}${F} &> /dev/null;then
       echo preloading ${F}...  
       scp -c -r ${F} ${LOCAL}
    fi
    TARGS=${LOCAL}`python -c "print '${F}'.split('/')[-1]"` ${TARGS}
done
SRC=`python -c "print '${ARGS[0]}'.split('/')[:-1]"`



echo lst_pipeline.sh $TARGS
/data1/paper/arp/scripts/lst_pipeline.sh $TARGS
for F in $TARGS; do
    echo Copying ${F}ddrxda to $SRC
    scp ${F}ddrxda ${SRC} && rm -r ${F}ddrxda
done
