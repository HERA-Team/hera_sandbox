#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=12G
#$ -l paper
#$ -N CALIBRATE_128
ARGS=`pull_args.py $*`

### XX-POL ###
POL='xx'
FIRSTCAL_S1E1='/data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz' # xx-pol
FIRSTCAL_S1E2='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz' # xx-pol
FIRSTCAL_S2E2='/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz' # xx-pol
FIRSTCAL_S2E3='/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz' # xx-pol
FIRSTCAL_S2E4='/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz' # xx-pol
FIRSTCAL_S2E5='/data4/paper/2014EoR/Analysis/ProcessedData/epoch5/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz' # xx-pol
BADANTS_S1E1='2,10,15,22,31,33,42,43,47,58,64,72,91,97,105,107,8,34,84,100' # xx-pol
BADANTS_S1E2='8,34,84,100' # xx-pol
BADANTS_S2E2='8,16,24,34,38,53,63,74,85' # xx-pol
BADANTS_S2E3='8,13,24,26,34,37,38,85,107' # xx-pol
BADANTS_S2E4='3,4,5,7,8,13,15,16,17,24,28,29,32,34,36,37,38,40,48,51,52,55,56,60,62,68,69,74,76,77,82,84,85,92,93,94,100,101,103,107,109,110' # xx-pol
BADANTS_S2E5='' # xx-pol

### YY-POL ###
#POL='yy'
#FIRSTCAL_S1E1='/data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz' # yy-pol
#FIRSTCAL_S1E2='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz' # yy-pol
#FIRSTCAL_S2E2='/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz' # yy-pol
#FIRSTCAL_S2E3='/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz' # yy-pol
#FIRSTCAL_S2E4='/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz' # yy-pol
#FIRSTCAL_S2E5='/data4/paper/2014EoR/Analysis/ProcessedData/epoch5/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz' # yy-pol
#BADANTS_S1E1='2,10,15,22,31,33,42,43,47,58,64,72,91,97,105,107,100,56,7,84'
#BADANTS_S1E2='100,56,7,84' # yy-pol
#BADANTS_S2E2='7,34,56,81' # yy-pol
#BADANTS_S2E3='7,15,16,17,26,34,56,82,107' # yy-pol
#BADANTS_S2E4='3,4,5,6,7,8,13,15,16,17,19,24,26,28,29,30,32,34,35,36,37,38,40,44,48,51,52,54,55,56,60,62,68,69,74,76,77,82,84,85,86,87,91,92,93,94,100,101,103,107,109,110' # yy-pol
#BADANTS_S2E5='' # yy-pol

#############################################################

echo ${ARGS}

for f in ${ARGS}; do
    : '
    ### SEASON 1 ###
    if ((${f:33:7} < 2456679 )); then
        echo working on ${f}, which is in Epoch 1...
        # Firstcal
        ~/capo/omni/omni_apply.py -p ${POL} --firstcal --fcfile ${FIRSTCAL_S1E1} ${f}
        # Omni Run
        ~/capo/omni/omni_run.py -p ${POL} -C psa6622_v003 --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v4_xtalk/ --ba ${BADANTS_S1E1} ${f:29:28}F
        # Omni Apply
        ~/capo/omni/omni_apply.py -p ${POL} --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba ${BADANTS_S1E1} ${f:29:28}F 
        echo deleting ${f:29:28}F
        rm -r ${f:29:28}F
    fi
    if (( ${f:33:7} > 2456678 )); then
        echo working on ${f}, which is in Epoch 2...
        # Firstcal
        ~/capo/omni/omni_apply.py -p ${POL} --firstcal --fcfile ${FIRSTCAL_S1E2} ${f}
        # Omni Run
        ~/capo/omni/omni_run.py -p ${POL} -C psa6622_v003 --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/ --ba ${BADANTS_S1E2} ${f:29:28}F
        # Omni Apply
        ~/capo/omni/omni_apply.py -p ${POL} --omnipath /data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba ${BADANTS_S1E2} ${f:29:28}F
        echo deleting ${f:29:28}F
        rm -r ${f:29:28}F
    fi
    '
    ### SEASON 2 ###
    if (( ${f:26:7} > 2456836 && ${f:26:7} < 2456875 )); then
        echo working on ${f}, which is in Epoch 2...
        # Firstcal
        ~/capo/omni/omni_apply.py -p ${POL} --ba ${BADANTS_S2E2} --firstcal --fcfile ${FIRSTCAL_S2E2} ${f}
        # Omni Run
        ~/capo/omni/omni_run.py -p ${POL} -C psa6622_v003 --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/ --ba ${BADANTS_S2E2} ${f:34:28}F
        # Omni Apply
        ~/capo/omni/omni_apply.py -p ${POL} --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba ${BADANTS_S2E2} ${f:34:28}F
        echo deleting ${f:34:28}F
        rm -r ${f:34:28}F
    fi
    if (( ${f:26:7} > 2456881 && ${f:26:7} < 2456929 )); then
        echo working on ${f}, which is in Epoch 3...
        # Firstcal
        ~/capo/omni/omni_apply.py -p ${POL} --ba ${BADANTS_S2E3} --firstcal --fcfile ${FIRSTCAL_S2E3} ${f}
        # Omni Run
        ~/capo/omni/omni_run.py -p ${POL} -C psa6622_v003 --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/ --ba ${BADANTS_S2E3} ${f:34:28}F
        # Omni Apply
        ~/capo/omni/omni_apply.py -p ${POL} --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba ${BADANTS_S2E3} ${f:34:28}F
        echo deleting ${f:34:28}F
        rm -r ${f:34:28}F
    fi 
    if (( ${f:26:7} > 2456941 && ${f:26:7} < 2456995 )); then
        echo working on ${f}, which is in Epoch 4...
        # Firstcal
        ~/capo/omni/omni_apply.py -p ${POL} --ba ${BADANTS_S2E4} --firstcal --fcfile ${FIRSTCAL_S2E4} ${f}
        # Omni Run
        ~/capo/omni/omni_run.py -p ${POL} -C psa6622_v003 --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/ --ba ${BADANTS_S2E4} ${f:34:28}F
        # Omni Apply
        ~/capo/omni/omni_apply.py -p ${POL} --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba ${BADANTS_S2E4} ${f:34:28}F
        echo deleting ${f:34:28}F
        rm -r ${f:34:28}F
    fi    
if (( ${f:26:7} > 2456996 )); then
        echo working on ${f}, which is in Epoch 5...
        # Firstcal
        ~/capo/omni/omni_apply.py -p ${POL} --ba ${BADANTS_S2E5} --firstcal --fcfile ${FIRSTCAL_S2E5} ${f}
        # Omni Run
        ~/capo/omni/omni_run.py -p ${POL} -C psa6622_v003 --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/ --ba ${BADANTS_S2E5} ${f:34:28}F
        # Omni Apply
        ~/capo/omni/omni_apply.py -p ${POL} --omnipath /data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba ${BADANTS_S2E5} ${f:34:28}F
        echo deleting ${f:34:28}F
        rm -r ${f:34:28}F
    fi    
done
