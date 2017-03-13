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
    #if ((${f:33:7} < 2456679 )); then
    #    echo working on ${f}, which is in Epoch 1...
    #    ~/capo/omni/omni_apply.py --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v3_xtalk/%s.npz -p yy --xtalk --ubls="0,2;0,1" --ba 2,10,15,22,31,33,42,43,47,58,64,72,91,97,105,107,100,7,56,84 ${f} #S1E1yy
    #fi
    
    if (( ${f:33:7} > 2456678 )); then
        echo working on ${f}, which is in Epoch 2...
        #~/capo/omni/omni_apply.py --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/%s.npz -p yy --xtalk --ubls="0,2;0,1" --ba 100,7,56,84 ${f} #S1E2yy
        ~/capo/omni/omni_apply.py --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk_removedegenOFF/%s.npz -p xx --xtalk --ubls="0,2;0,1" --ba 34,84,100 ${f} #S1E2xx
    fi
done
