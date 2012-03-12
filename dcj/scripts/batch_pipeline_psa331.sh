#$ -S /bin/bash
#$ -N psa332_pipe
#$ -j y
#$ -o grid_output

ARGS=`pull_args.py $*`

for ARG in $ARGS; do
#   echo correct_psa332_v005.py $ARG
#   correct_psa332_v005.py $ARG
   
   echo apply_bp.py ${ARG}
   apply_bp.py ${ARG}
   
   echo ~/scripts/xrfi_simple.py -c 0_143,323,378_387,850_852,900_1023 -t 53 ${ARG}*b
   ~/scripts/xrfi_simple.py -c 0_143,323,378_387,850_852,900_1023 -t 53 ${ARG}*b
   
   echo ~/scripts/flag_ant.py -a 0_1,2_3,4_5,6_7,8_9,10_11,12_13,14_15,16_17,18_19,20_21,22_23,24_25,26_27,28_29,30_31 ${ARG}*R
   ~/scripts/flag_ant.py -a 0_1,2_3,4_5,6_7,8_9,10_11,12_13,14_15,16_17,18_19,20_21,22_23,24_25,26_27,28_29,30_31 ${ARG}*R
done
echo $TASK_ID Finished