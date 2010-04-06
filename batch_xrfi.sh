#$ -S /bin/bash
ARGS=`pull_args.py $*`
#xrfi.py -m val -n 2.5 -c 0_49,200_255 $ARGS
xrfi.py -m val -n 2 -c 0_239,720_1023 $ARGS
