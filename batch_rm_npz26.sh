#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N rm_npz
#$ -o grid_output/

ARGS=`pull_args.py $*`

for ARG in $ARGS ; do
    /data1/paper/arp/scripts/rm_npz26.py -C psa331_v009_gc $ARG -d 45 -r 45
done
