#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo phs2src.py -C pgb015_v004 -s crab $ARGS
phs2src.py -C pgb015_v004 -s cas $ARGS
