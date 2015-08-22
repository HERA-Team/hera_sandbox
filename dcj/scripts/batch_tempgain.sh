#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=4G
#$ -j y
#$ -l paper
ARGS=`~/scripts/pull_args.py $*`
echo ~/scripts/tempgain.py --gom=1 --H_balun=-0.04 --H_cable=-0.018 $ARGS
~/scripts/tempgain.py --gom=1  --H_balun=-0.04 --H_cable=-0.018 $ARGS
