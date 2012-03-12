#$ -S /bin/bash
ARGS=`pull_args.py $*`
echo tempgain.py --H_balun=-0.0158 --H_cable=-0.0044 $ARGS
tempgain.py --gom=1 --H_balun=-0.0158 --H_cable=-0.0044 $ARGS
