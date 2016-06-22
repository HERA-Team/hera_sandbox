#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=12G
#$ -l scratch=9G
#$ -l paper
#$ -N OMNI_RUN
ARGS=`pull_args.py $*`

echo ${ARGS}

for f in ${ARGS}; do
    echo working on ${f}...
    ~/capo/omni/omni_run.py --calpar /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/calpar_first_cal_epoch3yy.p --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2_xtalk/ --ba 0,3,5,7,12,13,14,15,16,17,23,24,26,28,29,32,34,36,38,40,46,50,52,54,55,56,57,59,60,62,68,69,72,74,76,77,79,84,85,92,93,100,103,107,110 -p yy -C psa6622_v003 ${f}
    #~/capo/omni/omni_run.py --calpar /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/calpar_first_cal_epoch3yy.p --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2_xtalk/ --ba 5,7,8,15,16,17,24,26,27,28,29,37,38,46,48,50,51,53,55,63,68,69,72,74,76,77,82,83,84,85,92,107,110 -p xx -C psa6622_v003 ${f}
done
