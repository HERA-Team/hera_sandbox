#$ -S /bin/bash
#$ -V
#$ -cwd

ARGS=`pull_args.py $*`
CALFILE=psa6240_v001
for f in $ARGS; do
    echo apply_cal.py -C $CALFILE $ARGS/*.uvcRREcAzx
    apply_cal.py -C $CALFILE $ARGS/*.uvcRREcAzx
done
