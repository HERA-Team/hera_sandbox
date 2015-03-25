#$ -S /bin/bash
#$ -j y
#$ -o grid_output
#$ -e grid_output
#$ -cwd
#$ -V
ARGS=`pull_args.py $*`
#echo $ARGS

for FILE in $ARGS; do
    echo $FILE
done
