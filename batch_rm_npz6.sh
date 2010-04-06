#$ -S /bin/bash
ARGS=`pull_args.py $*`

rm_npz6.py -C pgb015_v005 $ARGS
