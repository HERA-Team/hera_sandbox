#$ -S /bin/bash
#$ -V
#$ -cwd
DIRS=`pull_args.py $*`
N=84
for f in $DIRS; do
    echo adjust_uv_filesize.py -n $N $f/*uvcRREcA
    adjust_uv_filesize.py -n $N $f/*uvcRREcA
done

