#$ -S /bin/bash
echo $PATH
ARGS=`pull_args.py $*`

for FILE in $ARGS; do
    xtalk3.py -o $FILE
done
