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
    omni_run.py --calpar /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/calpar_first_cal_epoch2xx.p --redinfo /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/redundantinfo_first_cal_epoch2xx.bin --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v2_alllsts/ ${f}
done

