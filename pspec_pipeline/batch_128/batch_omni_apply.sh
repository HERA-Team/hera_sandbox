#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=12G
#$ -l scratch=9G
#$ -l paper
#$ -N OMNI_APPLY
#$ -t 1:100
source activate PAPER
ARGS=`pull_args.py $*`

echo ${ARGS}

for f in ${ARGS}; do
    echo working on ${f}...
    ~/ReposForCanopy/capo/omni/omni_apply.py --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/sak_workspace/%s.npz -p yy --xtalk ${f} #epoch2yy season 2 epoch 2
done
