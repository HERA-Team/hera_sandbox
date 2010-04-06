#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo correct_pgb015.py $ARGS
correct_pgb015.py -t ~/pgb015/temps $ARGS
