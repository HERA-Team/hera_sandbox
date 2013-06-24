#$ -S /bin/bash
#$ -V
#$ -cwd
ARGS=`pull_args.py $*`
echo red_cal_v003.py --calpos=0,0 $ARGS
#xtalk3.py $ARGS
#xARGS=`python -c "print ' '.join(map(lambda x: x+'x', '$ARGS'.split()))"`
#echo $xARGS
red_cal_v003.py --calpos=0,0 $ARGS
