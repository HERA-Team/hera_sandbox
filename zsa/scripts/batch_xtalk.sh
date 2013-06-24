#$ -S /bin/bash
#$ -V
#$ -cwd

ARGS=`pull_args.py $*`
echo xtalk3.py  $ARGS 
xtalk3.py  $ARGS 
