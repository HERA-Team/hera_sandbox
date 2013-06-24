#$ -S /bin/bash
#$ -V
#$ -cwd

ARGS=`pull_args.py $*`
CORRECT=correct_psa6240.py

echo $CORRECT $ARGS
$CORRECT $ARGS
