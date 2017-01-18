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
    if (( ${f:33:7} > 2456678 )); then
        echo working on ${f}, which is in E2...
        #~/capo/omni/omni_run.py  -p yy -C psa6622_v003 --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/ --fc2 /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk/*yy.uvcRREc.median.fc.npz --ba 100,7,56,84 ${f} #S1E2yy 
        ~/capo/omni/omni_run.py  -p xx -C psa6622_v003 --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk_removedegenOFF/ --fc2 /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v3_xtalk_removedegenOFF/*xx.uvcRREc.median.fc.npz --ba 34,84,100 ${f} #S1E2xx 
    fi 
    #if (( ${f:33:7} < 2456679 )); then
    #    echo working on ${f}, which is in E1...
    #    ~/capo/omni/omni_run.py  -p yy  -C psa6622_v003 --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v3_xtalk/ --fc2 /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v3_xtalk/*yy.uvcRREc.median.fc.npz --ba 2,10,15,22,31,33,42,43,47,58,64,72,91,97,105,107,100,7,56,84 ${f} #S1E1yy
    #fi
    #~/capo/omni/omni_run.py --calpar /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/calpar_first_cal_epoch1xx.p --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v2_xtalk/ --ba 2,10,14,15,16,22,26,27,28,31,33,34,38,42,43,46,47,50,53,58,64,72,74,84,91,97,105,107 -p xx -C psa6622_v003 ${f} #S1E1xx
    #~/capo/omni/omni_run.py --calpar /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/calpar_first_cal_epoch1yy.p --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v2_xtalk/ --ba 2,3,7,10,15,16,22,23,26,30,31,33,34,38,42,43,46,47,50,56,58,60,64,72,91,97,100,105,107 -p yy -C psa6622_v003 ${f} #S1E1yy
    #~/capo/omni/omni_run.py --calpar /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/calpar_first_cal_epoch3yy.p --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2_xtalk/ --ba 0,3,5,7,12,13,14,15,16,17,23,24,26,28,29,32,34,36,38,40,46,50,52,54,55,56,57,59,60,62,68,69,72,74,76,77,79,84,85,92,93,100,103,107,110 -p yy -C psa6622_v003 ${f} #S2E3yy
    #~/capo/omni/omni_run.py --calpar /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/calpar_first_cal_epoch3xx.p --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2_xtalk/ --ba 5,7,8,15,16,17,24,26,27,28,29,37,38,46,48,50,51,53,55,63,68,69,72,74,76,77,82,83,84,85,92,107,110 -p xx -C psa6622_v003 ${f} #S2E3xx
done
