#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -j y
#$ -l h_vmem=12G
#$ -l paper
#$ -N OMNI_RUN
ARGS=`pull_args.py $*`
BADANTS=`cat /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/2014_epoch2_bad_ants.txt`
CAL=psa6622_v002
echo ${ARGS}

for f in ${ARGS}; do
    echo working on ${f}...
    ~/src/capo/omni/omni_run.py \
        --calpar /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/calpar_first_cal_epoch2xx.p \
        --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v2.1/ --ba ${BADANTS} \
        -C ${CAL} ${f}
done

