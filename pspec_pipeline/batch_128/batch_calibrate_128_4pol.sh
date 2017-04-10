#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=30G
#$ -N CALIB_128_4POL

source activate PAPER

###
# Hand this script just a SINGLE POL of files
# Omnical in 4pol mode will find the others, 
# assuming they share the same directories
###
PATH2CAPO='/home/saulkohn/ReposForCanopy/capo'

ARGS=`pull_args.py $*`
POLS='xx,xy,yx,yy'
CAL='psa6622_v003'
CWD=$(pwd)

###
# each epoch requires a single median.fc.npz file for xx and yy
###

FCAL_s2e2_xx='/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz'
FCAL_s2e3_xx='/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz'
FCAL_s2e4_xx='/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/*xx.uvcRREc.median.fc.npz'

FCAL_s2e2_yy='/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz'
FCAL_s2e3_yy='/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz'
FCAL_s2e4_yy='/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/*yy.uvcRREc.median.fc.npz'

###
# join of bad antennae in xx and yy QA
###

BADANTS_s2e2='7,8,16,24,34,38,53,56,63,74,81,85' 
BADANTS_s2e3='7,8,16,24,34,38,56,85,107'
BADANTS_s2e4='7,8,16,17,24,27,28,29,34,38,50,54,68,85,103'

#######################
# Season 2 Processing #
#######################

###
# Firstcal run: create *F files for linear pols and simply copy for crosspols
###

echo 'Beginning firstcal application'

for f in ${ARGS}; do
    p=${f:52:2} #hardcoded for folio 2014EoR file structure
    fxx="${f/$p/xx}"
    fyy="${f/$p/yy}"
    fxy="${f/$p/xy}"
    fyx="${f/$p/yx}"
    
    echo cp -r ${fxy} ${CWD}/${fxy:34:28}F
    cp -r ${fxy} ${CWD}/${fxy:34:28}F
    
    echo cp -r ${fyx} ${CWD}/${fyx:34:28}F
    cp -r ${fyx} ${CWD}/${fyx:34:28}F
        
    if (( ${f:26:7} > 2456836 && ${f:26:7} < 2456875 )); then
        echo working on ${f}, which is in Epoch 2...
        echo ~/capo/omni/omni_apply.py -p xx --firstcal --fcfile ${FCAL_s2e2_xx} ${fxx}
        ${PATH2CAPO}/omni/omni_apply.py -p xx --firstcal --fcfile ${FCAL_s2e2_xx} ${fxx}

        echo ~/capo/omni/omni_apply.py -p yy --firstcal --fcfile ${FCAL_s2e2_yy} ${fyy}
        ${PATH2CAPO}/omni/omni_apply.py -p yy --firstcal --fcfile ${FCAL_s2e2_yy} ${fyy}
    fi
    
    if (( ${f:26:7} > 2456881 && ${f:26:7} < 2456929 )); then
        echo working on ${f}, which is in Epoch 3...
        
        echo ~/capo/omni/omni_apply.py -p xx --firstcal --fcfile ${FCAL_s2e3_xx} ${fxx}
        ${PATH2CAPO}/omni/omni_apply.py -p xx --firstcal --fcfile ${FCAL_s2e3_xx} ${fxx}

        echo ~/capo/omni/omni_apply.py -p yy --firstcal --fcfile ${FCAL_s2e3_yy} ${fyy}
        ${PATH2CAPO}/omni/omni_apply.py -p yy --firstcal --fcfile ${FCAL_s2e3_yy} ${fyy}
    fi

    if (( ${f:26:7} > 2456941 )); then
        echo working on ${f}, which is in Epoch 4...
        echo ~/capo/omni/omni_apply.py -p xx --firstcal --fcfile ${FCAL_s2e4_xx} ${fxx}
        ${PATH2CAPO}/omni/omni_apply.py -p xx --firstcal --fcfile ${FCAL_s2e4_xx} ${fxx}

        echo ~/capo/omni/omni_apply.py -p yy --firstcal --fcfile ${FCAL_s2e4_yy} ${fyy}
        ${PATH2CAPO}/omni/omni_apply.py -p yy --firstcal --fcfile ${FCAL_s2e4_yy} ${fyy}
    fi
done

echo 'FirstCalibration complete. Beginning Omnicalibration.'

###
# Firstcal processing complete. Now we can omnicalibrate.
# After omni_run and omni_apply, we remove the *F files,
# which are no longer required.
###

for f in ${ARGS}; do
    p=${f:52:2}
    fFc=${CWD}/${f:34:28}F
    if (( ${f:26:7} > 2456836 && ${f:26:7} < 2456875 )); then
        echo working on ${fFc}, which is in Epoch 2...
        echo ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/ --ba=${BADANTS_s2e2} ${fFc}
        ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/ --ba=${BADANTS_s2e2} ${fFc}
        echo ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba=${BADANTS_s2e2} ${fFc}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba=${BADANTS_s2e2} ${fFc}
    fi

    if (( ${f:26:7} > 2456881 && ${f:26:7} < 2456929 )); then
        echo working on ${fFc}, which is in Epoch 3...
        echo ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/ --ba=${BADANTS_s2e3} ${fFc}
        ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/ --ba=${BADANTS_s2e3} ${fFc}
        echo ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba=${BADANTS_s2e3} ${fFc}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba=${BADANTS_s2e3} ${fFc}
    fi
    
    if (( ${f:26:7} > 2456941 )); then
        echo working on ${fFc}, which is in Epoch 3...
        echo ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/ --ba=${BADANTS_s2e4} ${fFc}
        ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/ --ba=${BADANTS_s2e4} ${fFc}
        echo ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba=${BADANTS_s2e4} ${fFc}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/%s.npz --xtalk --ubls="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" --ba=${BADANTS_s2e4} ${fFc}
    fi
    
    echo Deleting F files for ${f}
    for pp in xx xy yx yy; do
        echo rm -r "${fFc/$p/$pp}"
        rm -r "${fFc/$p/$pp}"
    done
done
