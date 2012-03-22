#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo tempgain.py $ARGS
tempgain.py $ARGS
