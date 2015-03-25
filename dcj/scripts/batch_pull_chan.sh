#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l paper
#$ -l h_vmem=2G
#$ -j y
#$ -N slice
#$ -o grid_output
. /usr/global/paper/CanopyVirtualEnvs/shredddercanopyrc.sh
canopy-PAPER_Omni
FILES=`~/scripts/pull_args.py $*`
~/scripts/pull_chan.py $FILES
