#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=12G
#$ -l scratch=9G
#$ -l paper
#$ -N OMNI_RUN

source activate PAPER
ARGS=`pull_args.py $*`

echo ${ARGS}

for f in ${ARGS}; do
    echo working on ${f}...
    ~/capo/omni/omni_run.py --calpar=/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/calpar_first_cal_epoch2yy.p --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/sak_workspace/ --ba=3,7,13,15,16,23,26,34,38,46,50,51,72,81,107,110,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127  -p yy -C psa6622_v003 ${f}
done
