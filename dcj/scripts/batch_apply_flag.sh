#$ -S /bin/bash
#$ -N psa332_r
#$ -j y
#$ -o grid_output/
ARGS=`pull_args.py $*`
~/scripts/flag_chan.py -c 0_143,279_281,323_324,371_373,510,378_387,769,850_853,900_1023  $ARGS
