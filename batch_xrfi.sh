#$ -S /bin/bash
ARGS=`pull_args.py $*`
xrfi.py -m val -c 0_33,209_255 $ARGS
