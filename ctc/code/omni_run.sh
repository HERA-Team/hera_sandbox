#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=12G
#$ -l scratch=9G
#$ -l paper
ARGS=`pull_args.py $*`

for f in ${ARGS}; do
    echo working on ${f}...
    omni_run.py --calpar /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/calpar_first_cal_epoch2xx.p --redinfo /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/redundantinfo_first_cal_epoch2xx.bin --savepath /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/ ${f}
done

