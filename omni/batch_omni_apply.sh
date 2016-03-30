#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=12G
#$ -l scratch=9G
#$ -l paper
#$ -N OMNI_APPLY
ARGS=`pull_args.py $*`

echo ${ARGS}

for f in ${ARGS}; do
    echo working on ${f}...
    #~/capo/omni/omni_apply.py --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2_xtalk/%s.npz -p yy --xtalk ${f}
    #~/capo/omni/omni_apply.py --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/%s.npz -p xx --xtalk ${f} #epoch2xx season 1
    ~/capo/omni/omni_apply.py --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk_yy/%s.npz -p yy --xtalk ${f} #epoch2yy season 1
done
