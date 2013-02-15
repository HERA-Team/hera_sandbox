#$ -S /bin/bash
#$ -o grid_output
#$ -e grid_output
ARGS=`pull_args.py $*`

for FILE in $ARGS; do   
    data_redux.sh $FILE
done
