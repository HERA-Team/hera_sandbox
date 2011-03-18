#$ -S /bin/bash
#$ -N psa332_r
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd

CAL=psa455_v004_gc

ARGS=`pull_args.py $*`
bash_uv2ms.py -C ${CAL} ${ARGS}
