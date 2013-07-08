#$ -S /bin/bash
#$ -V
#$ -cwd

DIR=`pull_args.py $*`
CORRECT=correct_psa6240.py

echo $DIR

for dir in $DIR; do
    echo working on $dir
    echo $CORRECT $ARGS
    $CORRECT $dir/*.uvcRRE
done;
