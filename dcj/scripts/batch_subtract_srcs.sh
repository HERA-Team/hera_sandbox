#$ -S /bin/bash
#$ -N psa332_sub
#$ -j y
#$ -o grid_output

ARGS=`pull_args.py $*`
SRCS="hyd vir her sgr"

for ARG in $ARGS; do
    echo Processing $ARG, $SRCS
    uv_addsub.py --sub ${ARG} ${ARG}.e_hyd
    uv_addsub.py --sub ${ARG}d ${ARG}.e_vir 
    uv_addsub.py --sub ${ARG}dd ${ARG}.e_her 
    uv_addsub.py --sub ${ARG}ddd ${ARG}.e_sgr
done
ssh djm growlnotify -n PSA332_sub -d $JOBID -m "${SGE_JOB_ID}:${SGE_TASK_ID} is finished"