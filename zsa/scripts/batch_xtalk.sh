#$ -S /bin/bash
#$ -V
#$ -cwd

DIR=`pull_args.py $*`
for dir in $DIR; do 
    echo working on $dir 
    echo xtalk3.py  $dir/*.uvcRREcAz
    xtalk3.py  $dir/*.uvcRREcAz
done;
