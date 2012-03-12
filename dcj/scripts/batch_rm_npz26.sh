#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N rm_npz
#$ -o grid_output/

ARGS=`pull_args.py $*`

for ARG in $ARGS ; do
    /data1/paper/arp/scripts/rm_npz26.py -C psa455_v003_gc $ARG -d 15 -r 15
done
