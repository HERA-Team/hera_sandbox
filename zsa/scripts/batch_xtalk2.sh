#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper

DIR=`pull_args.py $*`
for dir in $DIR; do 
    echo working in $dir
    echo xtalk3.py $dir/*uvcRREcAz
    xtalk3.py $dir/*uvcRREcAz
done
    
 
