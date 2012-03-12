#$ -S /bin/bash
#$ -N psa332_F
#$ -j y
#$ -o grid_output/
#$ -V
#$ -cwd
ARGS=`pull_args.py $*`
#~/scripts/flag_ant.py -a 0_1,2_3,4_5,6_7,8_9,10_11,12_13,14_15,16_17,18_19,20_21,22_23,24_25,26_27,28_29,30_31, $ARGS
~/scripts/flag_ant.py -a `python -c "print ','.join(['%d_%d'%(i,i+9) for i in range(3,7)]+['%d_%d'%(i,i+10) for i in range(3,7)]+['%d_%d'%(i*2,i*2+1) for i in range(0,16)])"` $ARGS
#~/scripts/flag_wtf.py $ARGS