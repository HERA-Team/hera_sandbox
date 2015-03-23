#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper

ARGS=`pull_args.py $*`
for file in $ARGS; do 
    echo xtalk3.py $file
    xtalk3.py $file
done
 
