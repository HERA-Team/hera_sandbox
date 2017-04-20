#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=30G
#$ -N CALIB_128_4POL

source activate PAPER

############
# Preamble #
############

###
# Hand this script just a SINGLE POL of files
# Omnical in 4pol mode will find the others, 
# assuming they share the same directories
###
ARGS=`pull_args.py $*`

PATH2CAPO='/home/saulkohn/ReposForCanopy/capo'
OMNIPATH_s1e1='/data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v4_xtalk_4pol/'
OMNIPATH_s1e2='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk_4pol/'

OMNIPATH_s2e2='/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk_4pol_SAKtest/'
OMNIPATH_s2e3='/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_4polTEST/'
OMNIPATH_s2e4='/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_4polTEST/'

POLS='xx,xy,yx,yy'
CAL='psa6622_v003'
UBLS="0,3;0,2;0,1;-1,1;1,1;1,2;-1,2" #output baselines

CWD=$(pwd)

###
# each epoch requires a single median.fc.npz file for xx and yy
###

FCAL_s1e1_xx='/data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v4_xtalk/zen.2456638.17384.xx.uvcRREc.median.fc.npz'
FCAL_s1e2_xx='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/zen.2456680.17382.xx.uvcRREc.median.fc.npz'

FCAL_s2e2_xx='/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/zen.2456848.17386.xx.uvcRREc.median.fc.npz'
FCAL_s2e3_xx='/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/zen.2456913.17390.xx.uvcRREc.median.fc.npz'
FCAL_s2e4_xx='/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/zen.2456988.17389.xx.uvcRREc.median.fc.npz'

FCAL_s1e1_yy='/data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v4_xtalk/zen.2456638.17384.yy.uvcRREc.median.fc.npz'
FCAL_s1e2_yy='/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/zen.2456680.17382.yy.uvcRREc.median.fc.npz'

FCAL_s2e2_yy='/data4/paper/2014EoR/Analysis/ProcessedData/epoch2/omni_v4_xtalk/zen.2456848.17386.yy.uvcRREc.median.fc.npz'
FCAL_s2e3_yy='/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v4_xtalk/zen.2456913.17390.yy.uvcRREc.median.fc.npz'
FCAL_s2e4_yy='/data4/paper/2014EoR/Analysis/ProcessedData/epoch4/omni_v4_xtalk/zen.2456988.17389.yy.uvcRREc.median.fc.npz'

###
# join of bad antennae in xx and yy QA
###

BADANTS_s1e1='2,7,8,10,15,22,31,33,34,42,43,47,56,58,64,72,84,91,97,100,105,107'
BADANTS_s1e2='7,8,16,24,34,38,53,56,63,74,84,85,100'
BADANTS_s2e2='7,8,16,24,34,38,53,56,63,74,81,85' 
BADANTS_s2e3='7,8,16,24,34,38,56,85,107'
BADANTS_s2e4='7,8,16,17,24,27,28,29,34,38,50,54,68,85,103'

##############
# Processing #
##############

###
# Firstcal run: create *F files for linear pols and simply copy for crosspols
###

echo 'Beginning firstcal application'

for f in ${ARGS}; do
    ###
    # this is a pretty hacky method that makes the script easier to generalize the
    # rest of the script. It relies on no "z" being present in the rest of the
    # filepath.
    ###
    fQ=z`cut -d "z" -f 2 <<< "$f"`
    
    if ((${fQ:4:7} < 2456679 )); then
        echo working on ${f}, which is in Season 1 Epoch 1...
        echo ~/capo/omni/omni_apply.py -p ${POLS} --firstcal --fcfile=${FCAL_s1e1_xx},${FCAL_s1e1_yy} --ba=${BADANTS_s1e1} ${f}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --firstcal --fcfile=${FCAL_s1e1_xx},${FCAL_s1e1_yy} --ba=${BADANTS_s1e1} ${f}
    fi
    
    if ((${fQ:4:7} > 2456678 && ${fQ:4:7} < 2456836)); then
        echo working on ${f}, which is in Season 1 Epoch 2...
        echo ~/capo/omni/omni_apply.py -p ${POLS} --firstcal --fcfile=${FCAL_s1e2_xx},${FCAL_s1e2_yy} --ba=${BADANTS_s1e2} ${f}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --firstcal --fcfile=${FCAL_s1e2_xx},${FCAL_s1e2_yy} --ba=${BADANTS_s1e2} ${f}
    fi
    
    if (( ${fQ:4:7} > 2456836 && ${f:4:7} < 2456875 )); then
        echo working on ${f}, which is in Season 2 Epoch 2...
        echo ~/capo/omni/omni_apply.py -p ${POLS} --firstcal --fcfile=${FCAL_s2e2_xx},${FCAL_s2e2_yy} --ba=${BADANTS_s2e2} ${f}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --firstcal --fcfile=${FCAL_s2e2_xx},${FCAL_s2e2_yy} --ba=${BADANTS_s2e2} ${f}
    fi
    
    if (( ${fQ:4:7} > 2456881 && ${f:4:7} < 2456929 )); then
        echo working on ${f}, which is in Season 2 Epoch 3...
        echo ~/capo/omni/omni_apply.py -p ${POLS} --firstcal --fcfile ${FCAL_s2e3_xx},${FCAL_s2e3_yy} --ba=${BADANTS_s2e3} ${f}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --firstcal --fcfile ${FCAL_s2e3_xx},${FCAL_s2e3_yy} --ba=${BADANTS_s2e3} ${f}
    fi

    if (( ${fQ:4:7} > 2456941 )); then
        echo working on ${f}, which is in Season 2 Epoch 4...
        echo ~/capo/omni/omni_apply.py -p ${POLS} --firstcal --fcfile ${FCAL_s2e4_xx},${FCAL_s2e4_yy} --ba=${BADANTS_s2e4} ${f}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --firstcal --fcfile ${FCAL_s2e4_xx},${FCAL_s2e4_yy} --ba=${BADANTS_s2e4} ${f}
    fi
done

echo 'FirstCalibration complete. Beginning Omnicalibration.'

###
# Firstcal processing complete. Now we can omnicalibrate.
# After omni_run and omni_apply, we remove the *F files,
# which are no longer required.
###

for f in ${ARGS}; do
    fQ=z`cut -d "z" -f 2 <<< "$f"`
    p=${fQ:18:2}
    fFc=${CWD}/${fQ}F
    
    if ((${fQ:4:7} < 2456679 )); then
        echo working on ${fFc}, which is in Season 1 Epoch 1...
        echo ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s1e1} --ba=${BADANTS_s1e1} ${fFc}
        ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s1e1} --ba=${BADANTS_s1e1} ${fFc}
        echo ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s1e1}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s1e1} ${fFc}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s1e1}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s1e1} ${fFc}
    fi
    
    if ((${fQ:4:7} > 2456678 && ${fQ:4:7} < 2456836)); then
        echo working on ${fFc}, which is in Season 1 Epoch 2...
        echo ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s1e2} --ba=${BADANTS_s1e2} ${fFc}
        ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s1e2} --ba=${BADANTS_s1e2} ${fFc}
        echo ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s1e2}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s1e2} ${fFc}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s1e2}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s1e2} ${fFc}
    fi
    
    
    if (( ${fQ:4:7} > 2456836 && ${fQ:4:7} < 2456875 )); then
        echo working on ${fFc}, which is in Season 2 Epoch 2...
        echo ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s2e2} --ba=${BADANTS_s2e2} ${fFc}
        ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s2e2} --ba=${BADANTS_s2e2} ${fFc}
        echo ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s2e2}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s2e2} ${fFc}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s2e2}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s2e2} ${fFc}
    fi

    if (( ${fQ:4:7} > 2456881 && ${fQ:4:7} < 2456929 )); then
        echo working on ${fFc}, which is in Season 2 Epoch 3...
        echo ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s2e3} --ba=${BADANTS_s2e3} ${fFc}
        ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s2e3} --ba=${BADANTS_s2e3} ${fFc}
        echo ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s2e3}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s2e3} ${fFc}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s2e3}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s2e3} ${fFc}
    fi
    
    if (( ${fQ:4:7} > 2456941 )); then
        echo working on ${fFc}, which is in Season 2 Epoch 4...
        echo ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s2e4} --ba=${BADANTS_s2e4} ${fFc}
        ${PATH2CAPO}/omni/omni_run.py -p ${POLS} -C ${CAL} --minV --omnipath=${OMNIPATH_s2e4} --ba=${BADANTS_s2e4} ${fFc}
        echo ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s2e4}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s2e4} ${fFc}
        ${PATH2CAPO}/omni/omni_apply.py -p ${POLS} --omnipath=${OMNIPATH_s2e4}/%s.npz --xtalk --ubls=${UBLS} --ba=${BADANTS_s2e4} ${fFc}
    fi
    
    echo Deleting F files for ${f}
    for pp in xx xy yx yy; do
        echo rm -r "${fFc/$p/$pp}"
        rm -r "${fFc/$p/$pp}"
    done
done
