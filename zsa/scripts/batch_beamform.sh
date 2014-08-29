#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
DIRS=`pull_args.py $*`
CALFILE=psa6240_v003
#EXT=uvcRREcAzx
EXT=uvcRREcA
#minsep is to get rid of all north/south  baselines.
for dir in $DIRS; do
    echo beamform.py -p xx -C $CALFILE -s pic -a cross --minsep=95 $dir/*$EXT
    beamform.py -p xx -C $CALFILE -s pic -a cross --minsep=95 $dir/*$EXT
done
