#$ -S /bin/bash
#$ -V
#$ -cwd

ARGS=`pull_args.py $*`
CALFILE=psa6240_v000
echo apply_cal.py -C $CALFILE $ARGS
apply_cal.py -C $CALFILE $ARGS
